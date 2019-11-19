import requests
from requests.exceptions import HTTPError\

class gitlab():
    def __init__(self, url="https://gitlab.cern.ch", projectid="project", token="token"):
        self.url = url
        self.id = projectid
        self.id = self.id.replace("/","%2F")
        self.token = token
        self.header = {
                       'PRIVATE-TOKEN': token,
                       'user-agent':'python-request'
                      }
    def uploadFile(self, filename, returnContent=False):
        url = self.url+"/api/v4"+"/projects/"+self.id+"/uploads"
        print(url)
        r = requests.post(url, headers=self.header, files={'file':open(filename,'rb')}) 
        if returnContent:
            return r.content
        return r
    def postIssue(self, title, description=None, labels=None, returnContent=False):
        assert type(labels) in [list, str]
        url = self.url+"/api/v4"+"/projects/"+self.id+"/issues"
        attr = dict()
        attr["title"]=title
        if description:
            attr["description"]=description
        if type(labels) is list:
            attr["labels"]=",".join(labels)
        elif type(labels) is str:
            attr["labels"]=labels
        r = requests.post(url, headers=self.header, data=attr)
        if returnContent:
            return r.content
        return r

if __name__=="__main__":
    pass