pip install pandas
pip install synbiochem-py

export PYTHONPATH=$PYTHONPATH:.

python liv/ergothioneine/simulate.py \
	data/models/yeastGEM.xml \
	out/ergothioneine/ \
	ergothioneine \
	Swainston \
	Neil \
	"University of Liverpool" \
	neil.swainston@liverpool.ac.uk