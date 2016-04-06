from django.core.management.base import BaseCommand, CommandError
from django.db import models
from django.core.paginator import Paginator

def add_ids_to_compounds():
    """Assign project IDs to compounds by saving them"""
    CBHCompoundBatch = models.get_model("cbh_chembl_model_extension", "CBHCompoundBatch")
    Project = models.get_model("cbh_core_model", "Project")
    for p in list(Project.objects.all()):
        cs = CBHCompoundBatch.objects.filter(project=p).order_by("id")

        paginator = Paginator(cs, 100) # chunks of 1000

        for page in range(1, paginator.num_pages +1):
            bs = paginator.page(page).object_list
            for obj in bs:
                obj.save()
            
            # here you can do what you want with the row
            
            print "done page %d of %d" % (page ,paginator.num_pages)
        print "done %s" % p.name




class Command(BaseCommand):

    def handle(self, *args, **options):
        add_ids_to_compounds()







