"""episcaf cluster service: the real job broker behind the GUI tool.

`staging` holds the CPU staging logic shared with the local stdlib server; `app` is the FastAPI
service deployed on a cluster-access host and reached over a private SSH tunnel. See README.md.
"""
