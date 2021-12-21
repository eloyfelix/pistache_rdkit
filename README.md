# pistache_rdkit

Very dirty and simple C++ REST API example. Using defaults in most functions for now.

## To run it

### Build the Docker image and run a container

```bash
docker-compose up -d
```

### Query the service:

```bash
curl --request POST -d "O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2" http://localhost:9080/murckoScaffold
```

### To stop and delete the running container

```bash
docker-compose down
```
