from django.core.management.base import BaseCommand, CommandError
from django.db import models


def add_ids_to_compounds():
    """Assign project IDs to compounds by saving them"""
    CBHCompoundBatch = models.get_model("cbh_chembl_model_extension", "CBHCompoundBatch")
    Project = models.get_model("cbh_core_model", "Project")
    for p in Project.objects.all():
        cs = CBHCompoundBatch.objects.filter(project=p).order_by("id")
        for c in cs:
            c.save()
            print c.project_counter


class Command(BaseCommand):

    def handle(self, *args, **options):
        add_ids_to_compounds()







