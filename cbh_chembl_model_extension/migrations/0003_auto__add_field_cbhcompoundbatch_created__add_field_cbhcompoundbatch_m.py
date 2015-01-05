# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'CBHCompoundBatch.created'
        db.add_column(u'cbh_chembl_model_extension_cbhcompoundbatch', 'created',
                      self.gf('django.db.models.fields.DateTimeField')(default=datetime.datetime.now, blank=True),
                      keep_default=False)

        # Adding field 'CBHCompoundBatch.modified'
        db.add_column(u'cbh_chembl_model_extension_cbhcompoundbatch', 'modified',
                      self.gf('django.db.models.fields.DateTimeField')(default=datetime.datetime.now, blank=True),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'CBHCompoundBatch.created'
        db.delete_column(u'cbh_chembl_model_extension_cbhcompoundbatch', 'created')

        # Deleting field 'CBHCompoundBatch.modified'
        db.delete_column(u'cbh_chembl_model_extension_cbhcompoundbatch', 'modified')


    models = {
        u'cbh_chembl_model_extension.cbhcompoundbatch': {
            'Meta': {'object_name': 'CBHCompoundBatch'},
            'created': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now', 'blank': 'True'}),
            'ctab': ('django.db.models.fields.TextField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'custom_fields': (u'django_hstore.fields.DictionaryField', [], {}),
            'editable_by': (u'django_hstore.fields.DictionaryField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modified': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now', 'blank': 'True'}),
            'related_molregno_id': ('django.db.models.fields.IntegerField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'viewable_by': (u'django_hstore.fields.DictionaryField', [], {}),
            'warnings': (u'django_hstore.fields.DictionaryField', [], {})
        }
    }

    complete_apps = ['cbh_chembl_model_extension']