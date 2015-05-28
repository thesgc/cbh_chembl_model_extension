# -*- coding: utf-8 -*-
from django.db import models, connection
from django.contrib.contenttypes.models import ContentType
from django.contrib.auth.models import Permission, User, Group
import copy

from django_hstore import hstore

import random
from chembl_business_model.models import MoleculeDictionary
from filter_pains import detect_pains
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
from flowjs.models import FlowFile
from picklefield.fields import PickledObjectField
import re
from pybel import readstring
import json
from django.db.models.signals import post_save
from base64 import urlsafe_b64encode
import requests
from copy import deepcopy
from rdkit.Chem import  SDMolSupplier, AllChem, Draw, SanitizeMol, SanitizeFlags
from StringIO import StringIO
from django.template.defaultfilters import slugify
from cbh_chembl_model_extension.lookups import *


hstore.DictionaryField.register_lookup(KeyValuesAny)

hstore.DictionaryField.register_lookup(KeyValuesAll)

'''
Signal code - on a webauthed system, if a user is auto-created by a webauth login,
notify a superuser.
'''
if "django_webauth" in settings.INSTALLED_APPS:
    from django.db.models.signals import post_save
    from django.core.mail import send_mail

    def email_new_user(sender, instance, **kwargs):
        if kwargs["created"]:  # only for new users
            #need to email superusers to inform them a new user has logged in
            if instance.email:
                admin_users = find_superuser()
                email_from = 'no-reply-chemreg@chembiohub.ox.ac.uk'
                new_user_name = '%s %s' % (instance.first_name, instance.last_name)
                for admin in admin_users:
                    html_email = '<p>Hi %s, <br>A new user has logged onto the system via Webauth, %s. <br>You should add them to any applicable projects.<br>Thanks<br>ChemBio Hub ChemReg</p>' % (admin.first_name, new_user_name)
                    email_message = 'Hi %s, A new user has logged onto the system via Webauth, %s.You should add them to any applicable projects' % (admin.first_name, new_user_name)
                    send_mail('New Webauth User', email_message, email_from, [admin.email,], fail_silently=False, html_message=html_email)
                #we also need to email the user with a welcome message
                welcome_message = 'Welcome to Chemreg, %s! You have successfully logged in using Webauth and you will be added to projects by an admin user shortly.' % new_user_name
                html_welcome_message = '<p>Welcome to Chemreg, %s! <br>You have successfully logged in using Webauth and you will be added to projects by an admin user shortly.</p>' % new_user_name
                send_mail('Welcome to ChemReg', welcome_message, email_from, [instance.email,], fail_silently=False, html_message=html_welcome_message)



    def find_superuser():
        return User.objects.all().filter(is_active=True, is_superuser=True).exclude(email=None).exclude(email="")


    post_save.connect(email_new_user, sender=User)



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



def get_all_hstore_values(table,column, key, is_list=False, extra_where=" True"):
    '''Using an hstore query from the reference here
    http://www.youlikeprogramming.com/2012/06/mastering-the-postgresql-hstore-date-type/
    where project_id = {{project_id}}
    '''
    cursor = connection.cursor()
    sql = "SELECT DISTINCT {column} -> '{key}' FROM {table} where {column} -> '{key}' != '' and {extra_where};".format( **{"table":table, "key":key, "column": column,  "extra_where": extra_where})
    cursor.execute(sql)
    mytuple = cursor.fetchall()
    items = []
    data = [d[0] for d in mytuple]
    for d in data:
        if is_list:
            try:
                for elem in json.loads(d):
                    items.append(elem)
            except ValueError:
                items.append(d)
        else:
            items.append(d)
    return items

BONDS_WEDGED_SDF_PROP = '''>  <_drawingBondsWedged>
True

$$$$'''
class CBHCompoundBatchManager(hstore.HStoreManager):
    def blinded(self,project=None):
        '''Generate a batch with a blinded id'''
      
        return CBHCompoundBatch(project=project, blinded_batch_id="EMPTY_ID")


    def from_rd_mol(self, rd_mol, orig_ctab=None,smiles="", project=None, reDraw=None):
        '''Clean up the structures that come in from Smiles or from XLS or SDFs'''
        #Get a copy of the mol data
        moldata = rd_mol
        if orig_ctab is None:
            for name in moldata.GetPropNames():
                #delete the property names for the saved ctab
                moldata.ClearProp(name)
            ctab = Chem.MolToMolBlock(moldata)
            #make it into an sdf block
            ctab += "\n$$$$"
        else:
            not_first_lists = orig_ctab.split("\n")
            headerlines = True
            lines = []
            reallines =0
            for line in not_first_lists:
                if "V2000" in line:
                    headerlines = False
                if headerlines:
                    pass
                else:
                    lines.append(line)
                    reallines += 1
                if "END" in line:
                    break
            if reallines == 2:
                raise Exception("blank_molfile")
            all_lines =["","","",] + lines + [BONDS_WEDGED_SDF_PROP,]

            ctab = "\n".join(all_lines)
        batch = CBHCompoundBatch(ctab=ctab, original_smiles=smiles, )
        batch.project_id = project.id
        return batch

    def get_all_keys(self, where=True, secondwhere=True):
        cursor = connection.cursor()
        cursor.execute("SELECT key, count(*) FROM (SELECT (each(custom_fields)).key FROM cbh_chembl_model_extension_cbhcompoundbatch WHERE %s) AS stat WHERE %s GROUP BY key  ORDER BY count DESC, key" % (where,secondwhere))
        # cursor.execute("select distinct k from (select skeys(custom_fields) as k from cbh_chembl_model_extension_cbhcompoundbatch where %s) as dt" % where )
        return cursor.fetchall()

    

    def index_new_compounds(self):
        cursor = connection.cursor()
        cursor.execute("INSERT INTO compound_mols (molregno , ctab) SELECT c.molregno, mol_from_ctab(molfile::cstring) ctab FROM compound_structures c LEFT OUTER JOIN compound_mols ON c.molregno = compound_mols.molregno WHERE is_valid_ctab(molfile::cstring) AND compound_mols.molregno is null;")
        return True




PROJECT_PERMISSIONS = (("viewer","Can View"),
                        ("editor","Can edit or add batches"),
                        ( "admin", "Can assign permissions"))

class ProjectPermissionManager(models.Manager):

    def sync_all_permissions(self):
        for perm in self.all():
            perm.sync_permissions()

    def get_user_permission(self,project_id, user, codenames, perms=None):
        '''Check the given users' permissions against a list of codenames for a project id'''
        if not perms:
            perms = user.get_all_permissions()
        codes = ["%d.%s" % (project_id, codename) for codename in codenames]
        matched = list(perms.intersection(codes))
        if len(matched) > 0:
            return True
        return False

    



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
            group_or_user.permissions.add(self.get_instance_permission_by_codename(codename))
        if type(group_or_user) == User:
            group_or_user.user_permissions.add(self.get_instance_permission_by_codename(codename))

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
    custom_field_config = models.ForeignKey("cbh_chembl_model_extension.CustomFieldConfig", related_name="project",null=True, blank=True, default=None, )
    is_default = models.BooleanField(default=False)
    
    class Meta:
        get_latest_by = 'created'

    def __unicode__(self):
        return self.name

    @models.permalink
    def get_absolute_url(self):
        return {'post_slug': self.project_key}



def sync_permissions(sender, instance, created, **kwargs):
    '''After saving the project make sure it has entries in the permissions table'''
    if created is True:
        instance.sync_permissions()
        instance.make_editor(instance.created_by)

post_save.connect(sync_permissions, sender=Project, dispatch_uid="proj_perms")





class CustomFieldConfig(TimeStampedModel):
    name = models.CharField(unique=True, max_length=50)
    created_by = models.ForeignKey("auth.User")
    schemaform = models.TextField(default = "", null=True, blank=True, )
    def __unicode__(self):
        return self.name

    

class CBHCompoundMultipleBatch(TimeStampedModel):
    '''Holds a list of batches'''
    created_by = models.CharField(max_length=50, db_index=True, null=True, blank=True, default=None)
    project = models.ForeignKey(Project, null=True, blank=True, default=None)    
    uploaded_data = PickledObjectField()
    uploaded_file = models.OneToOneField(FlowFile, null=True, blank=True, default=None)
    saved = models.BooleanField(default=False)
    #batches = models.ForeignKey(CBHCompoundBatch, null=True, default=None)







class PinnedCustomField(TimeStampedModel):
    TEXT = "text"
    TEXTAREA = "textarea"
    UISELECT = "uiselect"
    INTEGER  = "integer"
    NUMBER = "number"
    UISELECTTAG  = "uiselecttag"
    UISELECTTAGS = "uiselecttags"
    CHECKBOXES = "checkboxes"
    PERCENTAGE = "percentage"
    DATE = "date"

    FIELD_TYPE_CHOICES = {
                            "char" : {"name" : "Short text field", "data": { "type": "string" }},

                            TEXT : {"name" : "Short text field", "data": { "type": "string" }},
                            TEXTAREA: {"name" :"Full text", "data": { "type": "string" , "format" : "textarea"}},
                            UISELECT: {"name" :"Choice field", "data": { "type": "string" , "format" : "uiselect"}},
                            INTEGER: {"name" :"Integer field", "data": { "type": "integer"}},
                            NUMBER: {"name" :"Decimal field", "data": { "type": "number"}},
                            UISELECTTAG: {"name" : "Choice allowing create", "data":  { "type": "string", "format" : "uiselect"}},
                            UISELECTTAGS: {"name" : "Tags field allowing create" , "data": { "type": "array", "format" : "uiselect", "options": {
                                      "tagging": "tagFunction" ,
                                      "taggingLabel": "(adding new)",
                                      "taggingTokens": "",
                                 }}},
                            PERCENTAGE: {"name" :"Percentage field", "data": { "type": "number", "maximum" : 100.0, "minimum": 0.0}},
                            DATE:  {"name": "Date Field - today or past" , "data":{"type": "string",   "format": "date"}},
                        
                        }
    #                            CHECKBOXES : {"name": "Checkbox Fields", "data": {"type" : "array", "format": "checkboxes"}}


    field_key = models.CharField(max_length=50,  default="")
    name = models.CharField(max_length=50)
    description = models.CharField(max_length=1024, blank=True, null=True, default="")
    custom_field_config = models.ForeignKey("cbh_chembl_model_extension.CustomFieldConfig", related_name='pinned_custom_field')
    required = models.BooleanField(default=False)
    part_of_blinded_key = models.BooleanField(default=False, verbose_name="blind key")
    field_type = models.CharField(default="char", choices=((name, value["name"]) for name, value in FIELD_TYPE_CHOICES.items()), max_length=15, )
    allowed_values = models.CharField(max_length=1024, blank=True, null=True, default="")
    position = models.PositiveSmallIntegerField()

    def get_dropdown_list(self, projectKey):
        is_array = False
        if self.FIELD_TYPE_CHOICES[self.field_type]["data"]["type"] == "array":
            is_array=True
        db_items = get_all_hstore_values("cbh_chembl_model_extension_cbhcompoundbatch  inner join cbh_chembl_model_extension_project on cbh_chembl_model_extension_project.id = cbh_chembl_model_extension_cbhcompoundbatch.project_id ", 
            "custom_fields", 
            self.name, 
            is_list=is_array, 
            extra_where="cbh_chembl_model_extension_project.project_key ='%s'" % projectKey)
        return  [item for item in db_items]

    def get_allowed_items(self,projectKey):
        items = [item.strip() for item in self.allowed_values.split(",") if item.strip()]
        setitems = sorted(list(set(items + self.get_dropdown_list(projectKey))))
        testdata = [{"label" : item.strip(), "value": item.strip()} for item in setitems if item] 
        searchdata = [{"label" : "[%s] %s" % (self.name ,item.strip()), "value" : "%s|%s" % (self.name ,item.strip())} for item in setitems if item] 
        return (testdata, searchdata)

    class Meta:
        ordering = ['position']
        get_latest_by = 'created'






class CBHCompoundBatch(TimeStampedModel):
    '''Holds the batch information for an uploaded compound before it is fully registered'''
    ctab = models.TextField(null=True, blank=True, default=None)
    std_ctab = models.TextField(null=True, blank=True, default=None)
    canonical_smiles = models.TextField(null=True, blank=True, default=None)
    original_smiles = models.TextField(null=True, blank=True, default=None)
    editable_by =  hstore.DictionaryField() 
    uncurated_fields =  hstore.DictionaryField() 
    created_by = models.CharField(max_length=50, db_index=True, null=True, blank=True, default=None)
    #related_molregno_id = models.IntegerField(db_index=True,  null=True, blank=True, default=None)
    related_molregno = models.ForeignKey(MoleculeDictionary, null=True, blank=True, default=None, to_field="molregno", )
    standard_inchi = models.TextField(null=True, blank=True, default=None)
    standard_inchi_key = models.CharField(max_length=50,  null=True, blank=True, default=None)
    warnings =  hstore.DictionaryField() 
    properties = hstore.DictionaryField() 
    custom_fields =  hstore.DictionaryField() 
    errors = hstore.DictionaryField()
    multiple_batch_id = models.IntegerField(default=0)
    #multiple_batch_id = models.ForeignKey(CBHCompoundMultipleBatch, null=True, blank=True, default=None, to_field="id")
    objects = CBHCompoundBatchManager()
    project = models.ForeignKey(Project, null=True, blank=True, default=None)
    blinded_batch_id = models.CharField(default="",null=True, blank=True, max_length=12 )
    batch_number = models.IntegerField(default=-1,null=True, blank=True,)


    def get_uk(self,):
        return "%s__%d__%s" % (self.standard_inchi_key, self.project_id, "MOL")

    def save(self, *args, **kwargs):
        val = kwargs.pop("validate", True)

        # self.validate()
        for key in self.errors:
            raise ValidationError(key)
        if self.multiple_batch_id == 0:
            mb = CBHCompoundMultipleBatch.objects.create();
            self.multiple_batch_id = mb.id

        
        super(CBHCompoundBatch, self).save(*args, **kwargs)

    def validate(self, temp_props=True):         
        #self.set_pains_matches()
        self.standardise()
        # if temp_props:
        #     self.generate_temp_properties()
        

    # def set_pains_matches(self):

    #     self.warnings = detect_pains(self.std_ctab, {})

    def standardise(self):
        if self.canonical_smiles:
            return
        if not self.std_ctab:
            self.std_ctab = self.ctab
        if not self.standard_inchi:
            self.standard_inchi = inchiFromPipe(self.std_ctab, settings.INCHI_BINARIES_LOCATION['1.02'])
        if not self.standard_inchi:
            raise Exception("inchi_error")
        else:
            pass
            
            self.standard_inchi_key = InchiToInchiKey(self.standard_inchi)
                        


    def generate_structure_and_dictionary(self,chirality="1"):
        if self.id:
            #currently we dont update existing compound records
            pass
        else:
            if self.blinded_batch_id:
                
                uox_id = generate_uox_id()
                self.blinded_batch_id = uox_id
                self.save(validate=False)
                uox_id_lookup = ChemblIdLookup.objects.create(chembl_id=uox_id, 
                    entity_type="DOCUMENT",
                    entity_id=self.id )
            else:
                if not self.canonical_smiles :
                    pybelmol = readstring("inchi", self.standard_inchi)
                    self.canonical_smiles = pybelmol.write("can").split("\t")[0]
                    self.properties["cdxml"] = pybelmol.write("cdxml")

                    mol = MolFromInchi(self.standard_inchi)
                    if mol:
                        self.std_ctab = MolToMolBlock(mol, includeStereo=True)
                    inchi_key = self.standard_inchi_key
                    inchi = self.standard_inchi
                    if not self.related_molregno_id:
                        try:
                            moldict = MoleculeDictionary.objects.get(project=self.project, 
                                                                                structure_type="MOL",
                                                                                #chirality=chirality,
                                                                                structure_key=self.standard_inchi_key)
                        except ObjectDoesNotExist:

                            uox_id = generate_uox_id()
                            rnd = random.randint(-1000000000, -2 )
                            uox_id_lookup = ChemblIdLookup.objects.create(chembl_id=uox_id, entity_type="COMPOUND", entity_id=rnd)

                            moldict = MoleculeDictionary.objects.get_or_create(chembl=uox_id_lookup, 
                                                                                project=self.project, 
                                                                                structure_type="MOL",
                                                                                #chirality=chirality,
                                                                                structure_key=self.standard_inchi_key)[0]
                            uox_id_lookup.entity_id = moldict.molregno
                            uox_id_lookup.save()
                            structure = CompoundStructures(molecule=moldict,molfile=self.std_ctab, standard_inchi_key=inchi_key, standard_inchi=inchi)
                            structure.save()
                            generateCompoundPropertiesTask(structure)
                        self.related_molregno = moldict
                    self.save(validate=False)
