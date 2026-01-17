#!/bin/bash
# GPU Power Limit Script - Run at startup
# Sets RTX 3090 power limit to 200W for quieter operation

nvidia-smi -pm 1
nvidia-smi -pl 200
