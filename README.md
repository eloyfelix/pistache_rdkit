# pistache_rdkit

Very dirty and simple C++ REST API example. Using defaults in most functions for now.

## To run it

### Build the Docker image and run a container

```bash
docker-compose up -d
```

### Query the service:

```bash
curl --request POST -d "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1" http://localhost:9080/painsFilters
```

### To stop and delete the running container

```bash
docker-compose down
```
