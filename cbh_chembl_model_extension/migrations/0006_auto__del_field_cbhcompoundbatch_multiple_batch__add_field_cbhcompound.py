# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting field 'CBHCompoundBatch.multiple_batch'
        db.delete_column(u'cbh_chembl_model_extension_cbhcompoundbatch', 'multiple_batch_id')

        # Adding field 'CBHCompoundBatch.multiple_batch_id'
        db.add_column(u'cbh_chembl_model_extension_cbhcompoundbatch', 'multiple_batch_id',
                      self.gf('django.db.models.fields.IntegerField')(default=-1),
                      keep_default=False)


    def backwards(self, orm):
        # Adding field 'CBHCompoundBatch.multiple_batch'
        db.add_column(u'cbh_chembl_model_extension_cbhcompoundbatch', 'multiple_batch',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['cbh_chembl_model_extension.CBHCompoundMultipleBatch'], null=True, blank=True),
                      keep_default=False)

        # Deleting field 'CBHCompoundBatch.multiple_batch_id'
        db.delete_column(u'cbh_chembl_model_extension_cbhcompoundbatch', 'multiple_batch_id')


    models = {
        u'cbh_chembl_model_extension.cbhcompoundbatch': {
            'Meta': {'object_name': 'CBHCompoundBatch'},
            'canonical_smiles': ('django.db.models.fields.TextField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'created': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.CharField', [], {'default': 'None', 'max_length': '50', 'null': 'True', 'db_index': 'True', 'blank': 'True'}),
            'ctab': ('django.db.models.fields.TextField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'custom_fields': (u'django_hstore.fields.DictionaryField', [], {}),
            'editable_by': (u'django_hstore.fields.DictionaryField', [], {}),
            'errors': (u'django_hstore.fields.DictionaryField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modified': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now', 'blank': 'True'}),
            'multiple_batch_id': ('django.db.models.fields.IntegerField', [], {}),
            'original_smiles': ('django.db.models.fields.TextField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'properties': (u'django_hstore.fields.DictionaryField', [], {}),
            'related_molregno_id': ('django.db.models.fields.IntegerField', [], {'default': 'None', 'null': 'True', 'db_index': 'True', 'blank': 'True'}),
            'standard_inchi': ('django.db.models.fields.TextField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'standard_inchi_key': ('django.db.models.fields.CharField', [], {'default': 'None', 'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'std_ctab': ('django.db.models.fields.TextField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'viewable_by': (u'django_hstore.fields.DictionaryField', [], {}),
            'warnings': (u'django_hstore.fields.DictionaryField', [], {})
        },
        u'cbh_chembl_model_extension.cbhcompoundmultiplebatch': {
            'Meta': {'ordering': "('-modified', '-created')", 'object_name': 'CBHCompoundMultipleBatch'},
            'created': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.CharField', [], {'default': 'None', 'max_length': '50', 'null': 'True', 'db_index': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modified': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now', 'blank': 'True'}),
            'uploaded_data': ('picklefield.fields.PickledObjectField', [], {}),
            'uploaded_file': ('django.db.models.fields.related.OneToOneField', [], {'default': 'None', 'to': u"orm['flowjs.FlowFile']", 'unique': 'True', 'null': 'True', 'blank': 'True'})
        },
        u'flowjs.flowfile': {
            'Meta': {'object_name': 'FlowFile'},
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'identifier': ('django.db.models.fields.SlugField', [], {'unique': 'True', 'max_length': '255'}),
            'original_filename': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'state': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'total_chunks': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'total_chunks_uploaded': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'total_size': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'updated': ('django.db.models.fields.DateField', [], {'auto_now': 'True', 'blank': 'True'})
        }
    }

    complete_apps = ['cbh_chembl_model_extension']