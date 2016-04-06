from django.core.management.base import BaseCommand, CommandError
from django.db import models
from django.core.paginator import Paginator

import gc

def queryset_iterator(queryset, chunksize=1000):
    '''''
    Iterate over a Django Queryset ordered by the primary key

    This method loads a maximum of chunksize (default: 1000) rows in it's
    memory at the same time while django normally would load all rows in it's
    memory. Using the iterator() method only causes it to not preload all the
    classes.

    Note that the implementation of the iterator does not support ordered query sets.
    '''
    pk = 0
    last_pk = queryset.order_by('-pk')[0].pk
    queryset = queryset.order_by('pk')
    while pk < last_pk:
        for row in queryset.filter(pk__gt=pk)[:chunksize]:
            pk = row.pk
            yield row
        gc.collect()


def add_ids_to_compounds():
    """Assign project IDs to compounds by saving them"""
    CBHCompoundBatch = models.get_model("cbh_chembl_model_extension", "CBHCompoundBatch")
    Project = models.get_model("cbh_core_model", "Project")
    for p in list(Project.objects.all()):
        cs = CBHCompoundBatch.objects.filter(project=p).exclude(project_counter=-1).order_by("id")

        
        for obj in queryset_iterator(cs):

            obj.save()
            
            
        print "done %s" % p.name




class Command(BaseCommand):

    def handle(self, *args, **options):
        add_ids_to_compounds()







