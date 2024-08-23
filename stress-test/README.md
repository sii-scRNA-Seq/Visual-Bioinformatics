# Stress Testing

Still not clear if we can validate  many expected responses to a single emission with artillery yet:
https://github.com/artilleryio/artillery/discussions/1750
https://github.com/artilleryio/artillery/issues/389
https://github.com/artilleryio/artillery/issues/170

## Install Artillery
```
$ npm install -g artillery@2.0.20
```

## Running

### Local
```
$ artillery run --output stress-test/local-stress.json --environment local stress-test/stress-test.yaml
```

### Production
```
$ artillery run --output stress-test/prod-stress.json --environment production stress-test/stress-test.yaml
```

### HTML Reports
```
$ artillery report --output stress-test/[local|prod]-report.html stress-test/[local|prod]-stress.json
```
