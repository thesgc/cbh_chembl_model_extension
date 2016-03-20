# -*- coding: utf-8 -*-
from django.db import models, connection
from django.contrib.auth.models import User

from django_hstore import hstore
import json
import random
from chembl_business_model.models import MoleculeDictionary
from chembl_business_model.models import CompoundStructures
from rdkit import Chem
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem.rdmolfiles import MolToMolBlock
from django.core.exceptions import ValidationError
from django.conf import settings
import requests
from chembl_business_model.models import CompoundStructures
from chembl_business_model.models import MoleculeDictionary
from chembl_business_model.models import ChemblIdLookup
from chembl_business_model.tasks import generateCompoundPropertiesTask
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django_extensions.db.models import TimeStampedModel
import shortuuid
from chembl_business_model.utils import inchiFromPipe
from rdkit.Chem import InchiToInchiKey
from rdkit.Chem import MolFromInchi
from rdkit.Chem import MolToMolBlock
from flowjs.models import FlowFile
from picklefield.fields import PickledObjectField
from pybel import readstring
from django.db.models.signals import post_save
from cbh_chembl_model_extension.lookups import *
from copy import copy
from rdkit.Chem.AllChem import Compute2DCoords
import base64
import StringIO
from rdkit.Chem import  SDMolSupplier, MolToMolBlock, MolFromSmarts, SDMolSupplier, AllChem, Draw, SanitizeMol, SanitizeFlags,   AssignAtomChiralTagsFromStructure

hstore.DictionaryField.register_lookup(KeyValuesAny)

hstore.DictionaryField.register_lookup(KeyValuesAll)
hstore.DictionaryField.register_lookup(KeyValuesSingle)

'''
Signal code - on a webauthed system, if a user is auto-created by a webauth login,
notify a superuser.
'''
if "django_webauth" in settings.INSTALLED_APPS:
    from django.db.models.signals import post_save
    from django.core.mail import send_mail

    def email_new_user(sender, instance, **kwargs):
        if kwargs["created"]:  # only for new users
            # need to email superusers to inform them a new user has logged in
            if instance.email:
                if instance.username != instance.email:
                    #This is not an invited user
                    admin_users = find_superuser()
                    email_from = 'no-reply-chemreg@chembiohub.ox.ac.uk'
                    new_user_name = '%s %s' % (
                        instance.first_name, instance.last_name)
                    for admin in admin_users:
                        html_email = '<p>Hi %s, <br>A new user has logged onto the system via Webauth, %s. <br>You should add them to any applicable projects.<br>Thanks<br>ChemBio Hub ChemReg</p>' % (
                            admin.first_name, new_user_name)
                        email_message = 'Hi %s, A new user has logged onto the system via Webauth, %s.You should add them to any applicable projects' % (
                            admin.first_name, new_user_name)
                        send_mail('New Webauth User', email_message, email_from, [
                                  admin.email, ], fail_silently=False, html_message=html_email)
                    # we also need to email the user with a welcome message
                    welcome_message = 'Welcome to ChemReg, %s! You have successfully logged in using Webauth and you will be added to projects by an admin user shortly.' % new_user_name
                    html_welcome_message = '<p>Welcome to ChemReg, %s! <br>You have successfully logged in using Webauth and you will be added to projects by an admin user shortly.</p>' % new_user_name
                    send_mail('Welcome to ChemReg', welcome_message, email_from, [
                              instance.email, ], fail_silently=False, html_message=html_welcome_message)

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

    uox_id = "%s%s%s%s" % (
        settings.ID_PREFIX, two_letter, two_number, three_letter)
    try:
        ChemblIdLookup.objects.get(chembl_id=uox_id)
        return generate_uox_id()
    except ObjectDoesNotExist:
        return uox_id


BONDS_WEDGED_SDF_PROP = '''>  <_drawingBondsWedged>
True

$$$$'''


#-----------------------------------------------------------------------------------------------------------------------

def _parseMolData(data):
    suppl = SDMolSupplier()

    suppl.SetData(str(data), sanitize=False)
    data = [x for x in suppl if x]
    for x in data:
        if not x.HasProp("_drawingBondsWedged"):
            SanitizeMol(x)
        ctab = MolToMolBlock(x)
        ctablines = [item.split("0.0000") for item in ctab.split("\n") if "0.0000" in item]
        needs_redraw = 0
        for line in ctablines:
            if len(line) > 3:
                needs_redraw +=1
        if needs_redraw == len(ctablines):
             #check for overlapping molecules in the CTAB 
            SanitizeMol(x)
            Compute2DCoords(x)
            print "testr"
    return data

def _mols2imageStream(mols, f, format, size, legend, highlightMatch=None):
    highlights = None
    if highlightMatch:
        pattern = MolFromSmarts(highlightMatch)
        highlights = [mol.GetSubstructMatch(pattern) for mol in mols]
    kek = True
    if mols[0].HasProp("_drawingBondsWedged"):
        kek=False
    fit = False
    options = DrawingOptions() 
    if size >150:
        options.coordScale = 3
        options.bondLineWidth = 3.6
        options.dblBondOffset = 0.435
        options.atomLabelFontSize = 60
        fit = True
    
    image = Draw.MolsToGridImage(mols,molsPerRow=min(len(mols),4),subImgSize=(size*2,size*2),
                                     kekulize=kek,highlightAtomLists=highlights, fitImage=fit,
                                    options=options
    )
    image.save(f, format)

def _mols2imageString(mols,size,legend, format, recalc=False, highlightMatch=None):
    if not mols:
        return ''
 #   if recalc:
  #      _apply(mols, _computeCoords)
    imageData = StringIO.StringIO()
    for mol in mols:
        try:
            SanitizeMol(mol,sanitizeOps=SanitizeFlags.SANITIZE_ALL^SanitizeFlags.SANITIZE_CLEANUPCHIRALITY^Chem.SanitizeFlags.SANITIZE_SETCONJUGATION^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        except ValueError:
            return imageData.getvalue()
        AllChem.AssignAtomChiralTagsFromStructure(mol,replaceExistingTags=False)
    _mols2imageStream(mols, imageData, format, size, legend, highlightMatch=highlightMatch)
    return imageData.getvalue()

def _ctab2image(data,size,legend, recalc=True, highlightMatch=None):
    data = _mols2imageString(_parseMolData(data),size,legend, 'PNG', recalc=recalc, highlightMatch=highlightMatch)
    
    #if request.is_ajax:
    return base64.b64encode(data)


def set_images(batch):
    batch.bigimage = _ctab2image(copy(batch.ctab), 400, False, recalc=None)
    batch.image = _ctab2image(copy(batch.ctab),80,False, recalc=None)

class CBHCompoundBatchManager(hstore.HStoreManager):
    def get_image_for_assayreg(self, field, dpc, level):
        project = dpc.get("l0_permitted_projects")[0]
        if project[len(project) - 1] == "/":
            project = project[:-1]
        bits = project.split("/")
        id = bits[len(bits)-1]
        field_value = dpc[level]["project_data"].get(
            field["elasticsearch_fieldname"], None)
        #Fetch a batch with this UOx id if it is in the same project
        if field_value is not None:
            batches = CBHCompoundBatch.objects.filter(
                related_molregno__chembl__chembl_id=field_value,
                project_id=int(id))
            if batches.count() > 0:
                if batches[0].ctab:
                    return  _ctab2image(batches[0].ctab,75,False, recalc=None)
        return None

    def blinded(self, project=None):
        '''Generate a batch with a blinded id'''
        blinded_batch_id = "EMPTY_ID"
        return CBHCompoundBatch(project=project, blinded_batch_id=blinded_batch_id)

    def from_rd_mol(self, rd_mol, orig_ctab=None, smiles="", project=None, reDraw=None):
        '''Clean up the structures that come in from Smiles or from XLS or SDFs'''
        # Get a copy of the mol data
        moldata = rd_mol
        if orig_ctab is None and moldata:
            for name in moldata.GetPropNames():
                # delete the property names for the saved ctab
                moldata.ClearProp(name)
            ctab = Chem.MolToMolBlock(moldata)
            # make it into an sdf block
            ctab += "\n$$$$"
        else:
            not_first_lists = orig_ctab.split("\n")
            headerlines = True
            lines = []
            reallines = 0
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
            all_lines = ["", "", "", ] + lines + [BONDS_WEDGED_SDF_PROP, ]

            ctab = "\n".join(all_lines)

        batch = CBHCompoundBatch(ctab=ctab, original_smiles=smiles)
        set_images(batch)
        if project:
            batch.project_id = project.id
        return batch

    def get_all_keys(self, where=True, secondwhere=True):
        cursor = connection.cursor()
        cursor.execute(
            "SELECT key, count(*) FROM (SELECT (each(custom_fields)).key FROM cbh_chembl_model_extension_cbhcompoundbatch WHERE %s) AS stat WHERE %s GROUP BY key  ORDER BY count DESC, key" % (where, secondwhere))
        # cursor.execute("select distinct k from (select skeys(custom_fields) as k from cbh_chembl_model_extension_cbhcompoundbatch where %s) as dt" % where )
        return cursor.fetchall()

    def index_new_compounds(self):
        cursor = connection.cursor()
        cursor.execute(
            "INSERT INTO compound_mols (molregno , ctab) SELECT c.molregno, mol_from_ctab(molfile::cstring) ctab FROM compound_structures c LEFT OUTER JOIN compound_mols ON c.molregno = compound_mols.molregno WHERE is_valid_ctab(molfile::cstring) AND compound_mols.molregno is null;")
        return True


class CBHCompoundMultipleBatch(TimeStampedModel):

    '''Holds a list of batches'''
    created_by = models.CharField(
        max_length=50, db_index=True, null=True, blank=True, default=None)
    project = models.ForeignKey("cbh_core_model.Project", null=True, blank=True, default=None)
    uploaded_data = PickledObjectField()
    uploaded_file = models.ForeignKey(
        FlowFile, null=True, blank=True, default=None)
    saved = models.BooleanField(default=False)
    #batches = models.ForeignKey(CBHCompoundBatch, null=True, default=None)


class CBHCompoundBatch(TimeStampedModel):

    '''Holds the batch information for an uploaded compound before it is fully registered'''
    ctab = models.TextField(null=True, blank=True, default=None)
    std_ctab = models.TextField(null=True, blank=True, default=None)
    canonical_smiles = models.TextField(null=True, blank=True, default=None)
    original_smiles = models.TextField(null=True, blank=True, default=None)
    uncurated_fields = hstore.DictionaryField()
    image = models.TextField(default="")
    bigimage = models.TextField(default="")
    created_by = models.CharField(
        max_length=50, db_index=True, null=True, blank=True, default=None)
    created_by_id = models.IntegerField(null=True, blank=True, default=None)
    #related_molregno_id = models.IntegerField(db_index=True,  null=True, blank=True, default=None)
    related_molregno = models.ForeignKey(
        MoleculeDictionary, null=True, blank=True, default=None, to_field="molregno", )
    standard_inchi = models.TextField(null=True, blank=True, default=None)
    standard_inchi_key = models.CharField(
        max_length=50,  null=True, blank=True, default=None)
    warnings = hstore.DictionaryField()
    properties = hstore.DictionaryField()
    custom_fields = hstore.DictionaryField()
    multiple_batch_id = models.IntegerField(default=0)
    #multiple_batch_id = models.ForeignKey(CBHCompoundMultipleBatch, null=True, blank=True, default=None, to_field="id")
    objects = CBHCompoundBatchManager()
    project = models.ForeignKey("cbh_core_model.Project", null=True, blank=True, default=None)
    blinded_batch_id = models.CharField(
        default="", null=True, blank=True, max_length=12)
    batch_number = models.IntegerField(default=-1, null=True, blank=True,)




    def get_uk(self,):
        return "%s__%d__%s" % (self.standard_inchi_key, self.project_id, "MOL")

    def save(self, *args, **kwargs):
        val = kwargs.pop("validate", True)

        # self.validate()

        if self.multiple_batch_id == 0:
            mb = CBHCompoundMultipleBatch.objects.create()
            self.multiple_batch_id = mb.id
        super(CBHCompoundBatch, self).save(*args, **kwargs)
        if self.blinded_batch_id:
            uox_id_lookup = ChemblIdLookup.objects.get_or_create(chembl_id=self.blinded_batch_id,
                                                              entity_type="DOCUMENT",
                                                              entity_id=self.id)

        

    def validate(self, temp_props=True):
        self.standardise()

    def standardise(self):
        if self.canonical_smiles:
            return
        if not self.std_ctab:
            self.std_ctab = self.ctab
        if not self.standard_inchi:

            self.standard_inchi = inchiFromPipe(
                self.std_ctab, settings.INCHI_BINARIES_LOCATION['1.02'])
        if not self.standard_inchi:
            raise Exception("inchi_error")
        else:
            self.standard_inchi_key = InchiToInchiKey(
                self.standard_inchi.encode("ascii"))

    



def generate_structure_and_dictionary(batch):
    """
    Adding the structure data to a compound batch object
    """
    chirality="1"
    if batch.id:
        print "not updating"
        # currently we dont update existing compound records
    else:
        if batch.blinded_batch_id:

            uox_id = generate_uox_id()
            batch.blinded_batch_id = uox_id
            batch.save(validate=False)
            
        else:
            if not batch.canonical_smiles:
                try:
                    

                    pybelmol = readstring(
                        "mol", str(batch.ctab).encode("ascii"))
                    batch.canonical_smiles = pybelmol.write(
                        "can").split("\t")[0]
                    batch.properties["cdxml"] = pybelmol.write("cdxml")
                except:
                    pass
                try:
                    mol = MolFromInchi(
                        batch.standard_inchi.encode('ascii', 'ignore'))
                    if mol:
                        batch.std_ctab = MolToMolBlock(
                            mol, includeStereo=True)
                except:
                    pass
                inchi_key = batch.standard_inchi_key
                inchi = batch.standard_inchi
                if not batch.related_molregno_id:
                    try:
                        moldict = MoleculeDictionary.objects.get(project=batch.project,
                                                                 structure_type="MOL",
                                                                 # chirality=chirality,
                                                                 structure_key=batch.standard_inchi_key)
                    except ObjectDoesNotExist:

                        uox_id = generate_uox_id()
                        rnd = random.randint(-1000000000, -2)
                        uox_id_lookup = ChemblIdLookup.objects.create(
                            chembl_id=uox_id, entity_type="COMPOUND", entity_id=rnd)

                        moldict = MoleculeDictionary.objects.get_or_create(chembl=uox_id_lookup,
                                                                           project=batch.project,
                                                                           structure_type="MOL",
                                                                           # chirality=chirality,
                                                                           structure_key=batch.standard_inchi_key)[0]
                        uox_id_lookup.entity_id = moldict.molregno
                        uox_id_lookup.save()
                        structure = CompoundStructures(
                            molecule=moldict, molfile=batch.std_ctab, standard_inchi_key=inchi_key, standard_inchi=inchi)
                        structure.save()
                        if structure.molecule_id:
                            generateCompoundPropertiesTask(structure)
                    batch.related_molregno = moldict
                batch.save(validate=False)

    return batch



