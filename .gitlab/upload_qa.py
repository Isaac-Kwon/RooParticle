#
# Upload QA
# 

from gitlab import gitlab
from sys import argv
import json

print(argv)

g = gitlab(argv[1], argv[2], argv[3])
adr = g.uploadFile(argv[4],True)

adr = adr.decode('utf8').replace("'", '"')
adr = json.loads(adr)

content = """
Test QA

Produced MC File (ROOT): %s

"""%(adr["markdown"])

g.postIssue("Test QA", description=content, labels="Test QA, QA")




