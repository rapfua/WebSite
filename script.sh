#!/bin/bash

quarto render script.qmd 
shiny run --reload --launch-browser --port=0 /cloud/project/app.py