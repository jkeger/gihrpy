set -o xtrace

# Don't format the custom spacing etc in input object files
black . --exclude objects_
