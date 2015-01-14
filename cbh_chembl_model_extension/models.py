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
from rdkit.Chem.rdmolfiles import MolToMolBlock, MolFromMolBlock
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
from chembl_business_model.utils import inchiFromPipe
from rdkit.Chem import InchiToInchiKey
from rdkit.Chem import MolFromInchi
from rdkit.Chem import Kekulize
from rdkit.Chem import MolToMolBlock
from rdkit.Chem import MolFromMolBlock
from rdkit.Chem.rdmolfiles import MolToSmiles
from flowjs.models import FlowFile
from picklefield.fields import PickledObjectField



def generate_uox_id():
    two_letterg = shortuuid.ShortUUID()
    two_letterg.set_alphabet("ABCDEFGHJKLMNPQRSTUVWXYZ")
    two_letter = two_letterg.random(length=2)
    two_numberg = shortuuid.ShortUUID()
    two_numberg.set_alphabet("0123456789")
    two_number = two_numberg.random(length=2)
    three_letterg = shortuuid.ShortUUID()
    three_letterg.set_alphabet("ABCDEFGHJKLMNPQRSTUVWXYZ")
    three_letter = three_letterg.random(length=3)

    uox_id = "%s%s%s%s" % (settings.ID_PREFIX ,two_letter, two_number, three_letter )
    try:
        ChemblIdLookup.objects.get(chembl_id=uox_id)
        return generate_uox_id()
    except  ObjectDoesNotExist:
        return uox_id

class CBHCompoundBatchManager(hstore.HStoreManager):
    def from_rd_mol(self, rd_mol, smiles=""):
        batch = CBHCompoundBatch(ctab=Chem.MolToMolBlock(rd_mol), original_smiles=smiles)
        batch.validate(temp_props=False)
        return batch

    def get_all_keys(self, where=True):
        cursor = connection.cursor()
        cursor.execute("select distinct k from (select skeys(custom_fields) as k from cbh_chembl_model_extension_cbhcompoundbatch where %s) as dt" % where )
        return cursor.fetchall()


class CBHCompoundMultipleBatch(TimeStampedModel):
    '''Holds a list of batches'''
    created_by = models.CharField(max_length=50, db_index=True, null=True, blank=True, default=None)
    uploaded_data = PickledObjectField()
    uploaded_file = models.OneToOneField(FlowFile, null=True, blank=True, default=None)






class CBHCompoundBatch(TimeStampedModel):
    '''Holds the batch information for an uploaded compound before it is fully registered'''
    ctab = models.TextField(null=True, blank=True, default=None)
    std_ctab = models.TextField(null=True, blank=True, default=None)
    canonical_smiles = models.TextField(null=True, blank=True, default=None)
    original_smiles = models.TextField(null=True, blank=True, default=None)
    editable_by =  hstore.DictionaryField() 
    viewable_by =  hstore.DictionaryField() 
    created_by = models.CharField(max_length=50, db_index=True, null=True, blank=True, default=None)
    #related_molregno_id = models.IntegerField(db_index=True,  null=True, blank=True, default=None)
    related_molregno = models.ForeignKey(MoleculeDictionary, null=True, blank=True, default=None, to_field="molregno")
    standard_inchi = models.TextField(null=True, blank=True, default=None)
    standard_inchi_key = models.CharField(max_length=50,  null=True, blank=True, default=None)
    warnings =  hstore.DictionaryField() 
    properties = hstore.DictionaryField() 
    custom_fields =  hstore.DictionaryField() 
    errors = hstore.DictionaryField()
    multiple_batch_id = models.IntegerField(default=0)
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
        #managed=False


    def save(self, *args, **kwargs):
        val = kwargs.pop("validate", True)
        if val:
            self.validate()
            for key in self.errors:
                raise ValidationError(key)
        super(CBHCompoundBatch, self).save(*args, **kwargs)

    def validate(self, temp_props=True):
        self.warnings = {}
        self.errors = {}
             
        self.set_pains_matches()
        self.standardise()
        if temp_props:
            self.generate_temp_properties()
        

    def set_pains_matches(self):

        self.warnings = detect_pains(self.ctab, {})

    def standardise(self):
        warnings = []
        #testing without standardiser
        self.std_ctab = self.ctab
        # self.std_ctab = std(self.ctab,output_rules_applied=warnings, errors=self.errors)
        # for x, y in warnings:
        #     self.warnings[x] = y 
        self.standard_inchi = inchiFromPipe(self.std_ctab, settings.INCHI_BINARIES_LOCATION['1.02'])
        print self.standard_inchi
        mol = MolFromInchi(self.standard_inchi)
        self.canonical_smiles = MolToSmiles(mol)
        if not self.standard_inchi:
            self.errors["no_inchi"] = True
        self.standard_inchi_key = InchiToInchiKey(self.standard_inchi)
        self.std_ctab = MolToMolBlock(mol)
        self.warnings["hasChanged"] = self.original_smiles != self.canonical_smiles




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
        print ("generating structure")
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

