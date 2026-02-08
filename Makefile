# Colors for flavorrrr weeeeeeeeee~
RED		:= '\033[0;31m'
GREEN	:= '\033[0;32m'
YELLOW	:= '\033[1;33m'
BLUE	:= '\033[0;34m'
BOLD	:= '\033[1m'
NC		:= '\033[0m' # No Color

WASM_DIR	:= wasm
WASM_TARGET	:= $(WASM_DIR)/target/wasm32-unknown-unknown
WASM_OUT	:= website/wasm

PACKAGE_NAME	:= $(shell cd $(WASM_DIR) && cargo metadata --no-deps --format-version 1 | jq -r '.packages[0].name')
WASM_FILE	:= $(PACKAGE_NAME).wasm

.PHONY: help build watch deploy

help:
	@echo -e "$(BOLD)$(BLUE)build script$(NC)"
	@echo ""
	@echo -e "$(BOLD)COMMANDS:$(NC)"
	@echo -e "    $(GREEN)make build$(NC)    Build WASM (debug) and synth CDK"
	@echo -e "    $(GREEN)make watch$(NC)    Start a hot-reloading session"
	@echo    "                               served on http://localhost:8000/"
	@echo -e "    $(GREEN)make deploy$(NC)   Build everything (release) and deploy"
	@echo ""

build: cdk/node_modules # Build rust and synth CDK
	@echo -e "$(YELLOW)Building WASM...$(NC)"
	cd $(WASM_DIR) && cargo fmt
	cd $(WASM_DIR) && cargo build --target wasm32-unknown-unknown
	mkdir -p $(WASM_OUT)
	cp $(WASM_TARGET)/debug/$(WASM_FILE) $(WASM_OUT)/
	@echo -e "$(GREEN)Built WASM successfully, and copied $(WASM_FILE)"
	@echo -e "to $(BLUE)$(WASM_OUT)$(NC)\n"
	@echo -e "$(YELLOW)Building CDK...$(NC)"
	cd cdk && bunx cdk synth
	@echo -e "$(GREEN)Built CDK successfully!$(NC)\n"

deploy: cdk/node_modules # Build everything (release) and deploy
	@echo -e "$(YELLOW)Building WASM (release)...$(NC)"
	cd $(WASM_DIR) && cargo fmt
	cd $(WASM_DIR) && cargo build --target wasm32-unknown-unknown --release
	mkdir -p $(WASM_OUT)
	cp $(WASM_TARGET)/release/$(WASM_FILE) $(WASM_OUT)/
	@echo -e "$(GREEN)Built WASM successfully, and copied $(WASM_FILE)"
	@echo -e "to $(BLUE)$(WASM_OUT)$(NC)\n"
	@echo -e "$(YELLOW)Deploying...$(NC)"
	cd cdk && bunx cdk deploy --all
	@echo -e "$(GREEN)Deployed successfully!$(NC)\n"

cdk/node_modules: cdk/package.json
	cd cdk && bun install

CARGO_WATCH := $(WASM_DIR)/.tools/bin/cargo-watch

$(CARGO_WATCH):
	@echo -e "$(YELLOW)Fetching cargo-watch (this is a one-time process)$(NC)"
	mkdir -p -p $(WASM_DIR)/.tools
	cargo install --root $(WASM_DIR)/.tools --quiet cargo-watch

watch: $(CARGO_WATCH) # Start a hot-reloading session
	@echo -e "$(YELLOW)Building WASM...$(NC)"
	cd $(WASM_DIR) && cargo fmt
	cd $(WASM_DIR) && cargo build --target wasm32-unknown-unknown
	mkdir -p $(WASM_OUT)
	cp $(WASM_TARGET)/debug/$(WASM_FILE) $(WASM_OUT)/
	@echo -e "$(GREEN)Built WASM successfully, and copied $(WASM_FILE)"
	@echo -e "to $(BLUE)$(WASM_OUT)$(NC)\n"
	@echo -e "$(YELLOW)Starting watch mode...$(NC)"
	bunx serve --no-clipboard -p 8000 -L website &
	cd $(WASM_DIR) && $(CARGO_WATCH) -x "build --target wasm32-unknown-unknown" \
		-s "cp $(WASM_TARGET)/debug/$(WASM_FILE) $(WASM_OUT)/"
