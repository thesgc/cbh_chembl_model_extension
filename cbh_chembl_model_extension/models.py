# -*- coding: utf-8 -*-
from django.db import models, connection
from django_hstore import hstore
from chembl_business_model.models import MoleculeDictionary
from filter_pains import detect_pains
from standardiser.standardise import apply as std
from chembl_business_model.models import CompoundStructures
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Crippen
from rdkit.Chem.rdmolfiles import MolToMolBlock
from rdkit.Chem import rdMolDescriptors as Descriptors
from rdkit.Chem.SaltRemover import SaltRemover
from chembl_business_model.indigoWrapper import indigoObj
from django.core.exceptions import ValidationError
import StringIO
from django.db import IntegrityError
import requests
import base64
import decimal
from django.conf import settings
from django.db.models import Q
from chembl_business_model.utils import iterateModelRecords
from chembl_business_model.utils import iterateNModelRecords
from chembl_business_model.utils import ImageFromMolPP
from chembl_business_model.utils import cleanup
from chembl_business_model.models import CompoundImages
from chembl_business_model.models import CompoundProperties
from chembl_business_model.models import CompoundStructures
from chembl_business_model.models import CompoundRecords
from chembl_business_model.models import MoleculeHierarchy
from chembl_business_model.models import MoleculeDictionary
from chembl_business_model.models import ChemblIdLookup
from chembl_business_model.tasks import generateCompoundPropertiesTask

from chembl_business_model.models import Source
from django.db.models import Avg, Max, Min, Count
from django.db.models.fields import NOT_PROVIDED
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django_extensions.db.models import TimeStampedModel
import shortuuid

def generate_uox_id():
    uox_id = "%s_%s" % (settings.ID_PREFIX , shortuuid.ShortUUID().random(length=8))
    try:
        ChemblIdLookup.objects.get(chembl_id=uox_id)
        return generate_uox_id()
    except  ObjectDoesNotExist:
        return uox_id

class CBHBatchUpload(TimeStampedModel):
    '''file field upload for the individual file - one call to this service per file through FlowJS'''
    uploaded_file = models.FileField(upload_to="batchfiles")

class CBHCompoundBatchManager(hstore.HStoreManager):

    def get_all_keys(self, where=True):
        cursor = connection.cursor()
        cursor.execute("select distinct k from (select skeys(custom_fields) as k from cbh_chembl_model_extension_cbhcompoundbatch where %s) as dt" % where )
        return cursor.fetchall()



class CBHCompoundBatch(TimeStampedModel):
    '''Holds the batch information for an uploaded compound before it is fully registered'''
    ctab = models.TextField(null=True, blank=True, default=None)
    std_ctab = models.TextField(null=True, blank=True, default=None)
    editable_by =  hstore.DictionaryField() 
    viewable_by =  hstore.DictionaryField() 
    #related_molregno_id = models.IntegerField(null=True, blank=True, default=None)
    related_molregno = models.ForeignKey(MoleculeDictionary, null=True, blank=True, default=None, to_field="molregno")
    warnings =  hstore.DictionaryField() 
    properties = {} 
    custom_fields =  hstore.DictionaryField() 
    errors = {}
    #objects = hstore.HStoreManager()
    objects = CBHCompoundBatchManager()


    class Meta:
        '''In order to use as foreign key we set managed to false and set the migration to create the appropriate table
        The process is as follows - create the column that relates to chembl with _id as an integerfield
        Generate a set of migrations for that integer field
        Set managed = False
        run the migrations
        change the field to a ForeignKey - You now have a south migration for a model that relates to main chembl
        '''
        #pass
        managed=False


    def save(self, *args, **kwargs):
        self.validate()
        for key in self.errors:
            raise ValidationError(key)
        super(CBHCompoundBatch, self).save(*args, **kwargs)

    def validate(self):
        self.warnings = {}
        self.errors = {}
        try:
            
            self.set_pains_matches()
            self.standardise()
            self.generate_temp_properties()
        except:
            self.errors["invalid_molecule"] = True

    def set_pains_matches(self):

        self.warnings = detect_pains(self.ctab, {})

    def standardise(self):
        warnings = []
        self.std_ctab = std(self.ctab,output_rules_applied=warnings, errors=self.errors)
        for x, y in warnings:
            self.warnings[x] = y 


    def generate_temp_properties(self):
        prop = {}
        saltRemover = SaltRemover()
        mol = Chem.MolFromMolBlock(self.ctab)
        base = saltRemover.StripMol(mol)
        prop["hbd"] = Descriptors.CalcNumHBD(mol)
        prop["hba"] = Descriptors.CalcNumHBA(mol)
        prop["rtb"] = Descriptors.CalcNumRotatableBonds(mol)
        prop["alogp"] = Crippen.MolLogP(mol)
        prop["psa"] = Descriptors.CalcTPSA(mol)
        prop["full_mwt"] = Descriptors.CalcExactMolWt(mol)
        if base.GetNumAtoms():
            prop["mw_freebase"] = Descriptors.CalcExactMolWt(base)

        try:
            mol2 = indigoObj.loadMolecule(str(structure.molfile))
            prop["full_molformula"] = mol2.grossFormula()
        except:
            pass # TODO : handle this problem in smarter way

        self.properties = prop


    def generate_structure_and_dictionary(self):
        m = Chem.MolFromMolBlock(str(self.std_ctab))
        inchi = Chem.inchi.MolToInchi(m)
        inchi_key = Chem.inchi.InchiToInchiKey(inchi)
        try:
            structure = CompoundStructures.objects.get(standard_inchi_key=inchi_key)
            self.related_molregno_id = structure.molecule.molregno
            self.save()
        except ObjectDoesNotExist:
            uox_id_lookup = ChemblIdLookup.objects.create(chembl_id=generate_uox_id(), entity_type="COMPOUND")
            moldict = MoleculeDictionary.objects.get_or_create(chembl=uox_id_lookup)[0]
            uox_id_lookup.entity_id = moldict.molregno
            uox_id_lookup.save()
            structure = CompoundStructures(molecule=moldict,molfile=self.std_ctab, standard_inchi_key=inchi_key, standard_inchi=inchi)
            structure.save()
            self.related_molregno = moldict
            self.save()
        generateCompoundPropertiesTask(structure)

