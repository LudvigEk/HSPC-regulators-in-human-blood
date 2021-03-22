#!/bin/bash
python3 BCell_strategy_phase_II_CD34_to_NA.py &> phase2_no_CD34.log
python3 BCell_strategy_phase_II_for_missing_IgA.py &> phase2_no_IgA.log
python3 BCell_strategy_phase_II_upToQuadgate.py &> phase2_upToQuadgate.log
python3 BCell_strategy_phase_II_forplasmatransitonals.py &> phase2_forplasmatransitonals.log
