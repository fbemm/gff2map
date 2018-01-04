
import pprint

from BCBio.GFF import GFFExaminer

in_file = "your_file.gff"

examiner = GFFExaminer()

in_handle = open(in_file)

pprint.pprint(examiner.parent_child_map(in_handle))

in_handle.close()