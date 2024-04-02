#!/bin/zsh

# Source the zsh configuration files
if [ -f "$HOME/.zprofile" ]; then
    source "$HOME/.zprofile"
fi

if [ -f "$HOME/.zshrc" ]; then
    source "$HOME/.zshrc"
fi

./main

# Extras: 
    # source environment from zshrc (to run correctly)
    # Use cwd to set directory