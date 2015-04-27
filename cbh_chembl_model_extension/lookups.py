from django.db.models.lookups import In, Contains
from django_hstore.hstore import DictionaryField
from django_hstore.lookups import HStoreLookupMixin, get_cast_for_param






class KeyValues(HStoreLookupMixin,Contains):

    def get_statements(self, compiler, connection):
        lhs, lhs_params = self.process_lhs(compiler, connection)
        statements = []
        values = []
        for keyvalue in self.rhs:
            splitted = keyvalue.split("|")
            substs = {
                "field" : lhs,
                "key" : "\'%s\'" % splitted[0],
            }
            
            #Either match as part of an array or match exactly - check key is present so the exclude filter works
            statement = """((({field} -> {key}) like %s or ({field} -> {key}) = %s) and {field} ? {key})""".format(**substs) 
            statements.append(statement)
            values += ['%\\"' + splitted[1] + '\\"%', splitted[1]]

        return (statements, values)

class KeyValuesAll(KeyValues):
    lookup_name = 'kv_all'

    def as_postgresql(self, compiler, connection):
        statements, values = self.get_statements( compiler, connection)
        if len(statements) == 1:
            return (statements[0], values)
        else: 
            return (" and ".join( statements), values)

class KeyValuesAny(KeyValues):
    lookup_name = 'kv_any'

    def as_postgresql(self, compiler, connection):
        statements, values = self.get_statements( compiler, connection)
        if len(statements) == 1:
            return (statements[0], values)
        else: 
            return (" or ".join( statements), values)




