#!/bin/zsh

# Source the zsh configuration files
if [ -f "$HOME/.zprofile" ]; then
    source "$HOME/.zprofile"
fi

if [ -f "$HOME/.zshrc" ]; then
    source "$HOME/.zshrc"
fi

ctest -V