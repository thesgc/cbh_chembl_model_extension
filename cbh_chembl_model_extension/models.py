# -*- coding: utf-8 -*-
from django.db import models, connection
from django.contrib.contenttypes.models import ContentType
from django.contrib.auth.models import Permission, User, Group

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
import re
from pybel import readstring

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



PROJECT_PERMISSIONS = (("viewer","Can View"),
                        ("editor","Can edit or add batches"),
                        ( "admin", "Can assign permissions"))

class ProjectPermissionManager(models.Manager):

    def sync_all_permissions(self):
        for perm in self.all():
            perm.sync_permissions()

    def get_user_permission(self,project_id, user, codename):
        return user.has_perm("%d.%s" % (project_id, codename) )




class ProjectPermissionMixin(models.Model):
    '''The aim of this mixin is to create a permission content type and a permission model for a given project
    It allows for pruning the contnet types once the model is changed
    '''
    

    objects = ProjectPermissionManager()
    
    def get_project_key(self):
        return str(self.pk)

    def sync_permissions(self):
        '''first we delete the existing permissions that are not labelled in the model'''
        ct , created = ContentType.objects.get_or_create(app_label=self.get_project_key(), model=self, name=self.name)
        deleteable_permissions = Permission.objects.filter(content_type_id=ct.pk).exclude(codename__in=[perm[0] for perm in PROJECT_PERMISSIONS])
        deleteable_permissions.delete()
        for perm in PROJECT_PERMISSIONS:
            pm = Permission.objects.get_or_create(content_type_id=ct.id,codename=perm[0],name=perm[1])



    def get_contenttype_for_instance(self):
        ct = ContentType.objects.get(app_label=self.get_project_key(), model=self)
        return ct

    def delete_all_instance_permissions(self):
        '''for the pre delete signal'''
        deleteable_permissions = Permission.objects.filter(content_type_id=self.get_contenttype_for_instance().id)
        deleteable_permissions.delete()
                

    def get_instance_permission_by_codename(self, codename):
        pm = Permission.objects.get(codename=codename, content_type_id=self.get_contenttype_for_instance().id)
        return pm


    def _add_instance_permissions_to_user_or_group(self, group_or_user, codename):
        if type(group_or_user) == Group:
            group_or_user.group_permissions.add(self.get_instance_permission_by_codename(codename))
        if type(group_or_user) == User:
            group_or_user.user_permissions.add(self.get_instance_permission_by_codename(codename))
        group_or_user.save()

    def make_editor(self,group_or_user):
        self._add_instance_permissions_to_user_or_group(group_or_user, "editor")


    def make_viewer(self,group_or_user):
        self._add_instance_permissions_to_user_or_group(group_or_user, "viewer")


    def make_admin(self,group_or_user):
        self._add_instance_permissions_to_user_or_group(group_or_user, "admin")



    class Meta:
       
        abstract = True
        









class Project(TimeStampedModel, ProjectPermissionMixin):
    ''' Project is a holder for moleculedictionary objects and for batches'''
    name = models.CharField(max_length=50, db_index=True, null=True, blank=True, default=None)
    project_key = models.SlugField(max_length=50, db_index=True, null=True, blank=True, default=None, unique=True)
    created_by = models.ForeignKey("auth.User")

    class Meta:
        get_latest_by = 'created'

    def __unicode__(self):
        return self.name

    @models.permalink
    def get_absolute_url(self):
        return {'post_slug': self.project_key}



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
    project = models.ForeignKey(Project, null=True, blank=True, default=None)


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

        self.validate()
        for key in self.errors:
            raise ValidationError(key)
        #self.get_image_from_pipe()
        print self.properties
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
        pybelmol = readstring("inchi", self.standard_inchi)
        #pybel svg because the rdkit version does not support large organometallics
        self.properties["svg"] = pybelmol.write("svg")
        print self.properties["svg"]
        #pybel canonical smiles because the rdkit version does not support large organometallics
        self.canonical_smiles = pybelmol.write("can").split("\t")[0]
        #self.canonical_smiles = MolToSmiles(mol, canonical=True)
        
        if not self.standard_inchi:
            self.errors["no_inchi"] = True
        self.standard_inchi_key = InchiToInchiKey(self.standard_inchi)
        #rdkit molfile for rdkit database cartridge
        mol = MolFromInchi(self.standard_inchi)
        self.std_ctab = MolToMolBlock(mol)
        self.warnings["hasChanged"] = self.original_smiles != self.canonical_smiles




    # def get_image_from_pipe(self):
    #     '''Take a structure as a string ctab or string smiles and convert it to an svg
    #     format can be one of mol or smi 
    #     '''
    #     from subprocess import PIPE, Popen
    #     structure = self.original_smiles
    #     path = settings.OPEN_BABEL_EXECUTABLE
    #     p = Popen([path, "-ismi", "-xCe", "-osvg"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    #     a = p.communicate(input=str(structure))
    #     svg = a[0]

    #     self.properties["svg"] = re.sub(r'width="([0123456789\.]+)"\s+height="([0123456789.]+)"', 
    #         r'width="100%" viewbox="0 0 \1 \2" preserveAspectRatio="xMinYMin meet" version="1.1"', 
    #         svg)


    def generate_temp_properties(self):
        
        saltRemover = SaltRemover()
        mol = Chem.MolFromMolBlock(self.ctab)
        base = saltRemover.StripMol(mol)
        self.properties["hbd"] = Descriptors.CalcNumHBD(mol)
        self.properties["hba"] = Descriptors.CalcNumHBA(mol)
        self.properties["rtb"] = Descriptors.CalcNumRotatableBonds(mol)
        self.properties["alogp"] = Crippen.MolLogP(mol)
        self.properties["psa"] = Descriptors.CalcTPSA(mol)
        self.properties["full_mwt"] = Descriptors.CalcExactMolWt(mol)
        if base.GetNumAtoms():
            self.properties["mw_freebase"] = Descriptors.CalcExactMolWt(base)

        try:
            mol2 = indigoObj.loadMolecule(str(structure.molfile))
            self.properties["full_molformula"] = mol2.grossFormula()
        except:
            pass # TODO : handle this problem in smarter way


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

