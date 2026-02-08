# Welcome

A softbody simulation of a little blob

## Development


### Prerequisites

This project uses Typescript CDK so you're gonna need a JS runtime, like [Bun](https://bun.com/)

This project also uses Rust, so you're gonna need [that too](https://www.rust-lang.org/tools/install). 


### Project layout and wiring

This project has three main directories: `wasm/`, `website/`, and `cdk/`. Generally the entrypoint is `wasm/`.

#### `wasm/`
The `wasm/` directory holds all the rust code which gets compiled into WASM. The `build.sh` script uses the [`wasm-pack`](https://drager.github.io/wasm-pack/) utility to compile the rust and spit out the compiled stuff into `website/wasm/`.

#### `website/`
The `website/` directory has the website's files like `index.html` and `error.html`. `index.html` assumes that compiled WASM files are in `website/wasm/`, written there by `wasm-pack` as mentioned in the previous section.

#### `cdk/`
The `cdk/` directory has all of the infrastructure code, which is just 
- An S3 Bucket (holds all the website files like `.html` files, `.wasm` files, etc.)
- A [CloudFront Distribution](https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/distribution-working-with.html) which just serves up the files in da bucket
- Some DNS records injected into the HostedZone (which should be auto-created for you when you buy the domain) which point DNS to your CloudFront Distribution
- An SSL Certificate (makes da site HTTPS)

The `BucketDeployment` in `cdk/lib/website-stack.ts` takes files from `website/` and then plops them into the s3 bucket as-is. Dead simple. The BucketDeployment handles most of the complicated deployment stuff like CloudFront cache invalidation.

### How to build the project

1. Put the domain name in a `.env` file at the root of the project. so something like 
```text
  ├── cdk/
  ├── wasm/
  ├── website/
  ├── build.sh
    ...
  └── .env
```

Fill it with a variable `DOMAIN_NAME`
```
DOMAIN_NAME="your-website.com"
```
(This project assumes that you own this domain using [AWS Route53 domains](https://docs.aws.amazon.com/Route53/latest/DeveloperGuide/domain-register.html))

2. Build and deploy!
```bash
./build.sh --deploy
```

## Notes

If you wanna spin up the site locally:
```bash
./build.sh --watch
```

Then go to http://localhost:8000