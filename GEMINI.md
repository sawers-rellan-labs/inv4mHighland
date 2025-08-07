# Gemini CLI Instructions


This document outlines the operational guidelines and preferences for interacting with the Gemini CLI agent.

## General Interaction Principles

- **Tone:** Responses should be analytical and direct, avoiding overly enthusiastic or flattering language.
- **Code Generation:** When generating code, do not initially include report functions. Instead, suggest their inclusion and implement them only if explicitly requested.

## User's Preferred Languages

The user typically writes code in:
- R
- Bash
- Python

## R Coding Style Guidelines

Adhere to the following conventions when writing R code:

- **Indentation:** Use 2 spaces for indentation; tabs are not permitted.
- **Line Length:** Maintain a maximum of 80 characters per line.
- **Style:** Follow the tidyverse style guide, utilizing packages such as `dplyr`, `tidyr`, and `purrr`.
- **Function Export:** Export only functions explicitly marked with the `@export` roxygen tag.
- **Naming Conventions:** Use `snake_case` for all function and variable names.
- **Returns:** Prefer implicit returns (the last expression in a function is returned).
- **Input Validation:** Validate function inputs using specific validation functions.
- **Test Assertions:** Employ `expect_` functions for test assertions.
- **Documentation:** Document all exported functions with `roxygen2` comments, formatted using Markdown.
- **Parameter Type Validation:** Include parameter type validation at the beginning of functions.
- **Error Handling:** Follow tidyverse error handling patterns, providing informative messages.

folllow instructions from CLAUDE.md except thos regarding gemini-cli mcp 

