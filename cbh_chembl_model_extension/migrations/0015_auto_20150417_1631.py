# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


def migrate_multiple_batch_data(apps, schema_editor):
    # We can't import the Person model directly as it may be a newer
    # version than this migration expects. We use the historical version.
    Batch = apps.get_model("cbh_chembl_model_extension", "CBHCompoundBatch")

    for batch in Batch.objects.all():
        if batch.multiple_batch:
            batch.multiple_batch.created_by = batch.created_by
            batch.multiple_batch.project = batch.project
            batch.multiple_batch.save()



class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0014_auto_20150417_1629'),
    ]

    operations = [
        migrations.AddField(
            model_name='cbhcompoundmultiplebatch',
            name='project',
            field=models.ForeignKey(default=None, blank=True, to='cbh_chembl_model_extension.Project', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cbhcompoundmultiplebatch',
            name='saved',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
    ]
