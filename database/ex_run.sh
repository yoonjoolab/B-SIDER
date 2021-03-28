#!/bin/bash

# This script generates a complementary sequence database "example.db"
# and stores information of beta strands from 1chd.pdb

pdb=1chdA.pdb
pymol $pdb -cq B-SIDER_DB_builder.py -- $pdb example.db
