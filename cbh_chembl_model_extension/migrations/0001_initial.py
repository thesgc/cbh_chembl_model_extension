# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'CBHCompoundBatch'
        db.execute("create extension if not exists hstore")
        db.create_table(u'cbh_chembl_model_extension_cbhcompoundbatch', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('ctab', self.gf('django.db.models.fields.TextField')(default=None, null=True, blank=True)),
            ('editable_by', self.gf(u'django_hstore.fields.DictionaryField')()),
            ('viewable_by', self.gf(u'django_hstore.fields.DictionaryField')()),
            ('related_molregno_id', self.gf('django.db.models.fields.IntegerField')()),
            ('warnings', self.gf(u'django_hstore.fields.DictionaryField')()),
            ('custom_fields', self.gf(u'django_hstore.fields.DictionaryField')()),
        ))
        db.send_create_signal(u'cbh_chembl_model_extension', ['CBHCompoundBatch'])


    def backwards(self, orm):
        # Deleting model 'CBHCompoundBatch'
        db.delete_table(u'cbh_chembl_model_extension_cbhcompoundbatch')


    models = {
        u'cbh_chembl_model_extension.cbhcompoundbatch': {
            'Meta': {'object_name': 'CBHCompoundBatch'},
            'ctab': ('django.db.models.fields.TextField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'custom_fields': (u'django_hstore.fields.DictionaryField', [], {}),
            'editable_by': (u'django_hstore.fields.DictionaryField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'related_molregno_id': ('django.db.models.fields.IntegerField', [], {}),
            'warnings': (u'django_hstore.fields.DictionaryField', [], {}),
            'viewable_by': (u'django_hstore.fields.DictionaryField', [], {})
        }
    }

    complete_apps = ['cbh_chembl_model_extension']