from urllib.request import urlopen
import sys
import os
import time

fname = 'invivodec02' # file name
ckpoint = 'checkpoint.cpk'

delim = '|'  # file delimiter

# columns in file.  I don't know what the last one is
columns = ['tumor', 'strain', 'site', 'schedule', 'route', 'vehicle', 'parameter', 'NSC', 'unit', 'dose', 'T/C', 'endpoint', 'survived', 'x' ]


def read_dict(fname):
	""" read mapping dictionaries from the NCI data instructions """
	result = dict()
	with open(fname, 'r') as file:
		for line in file.readlines():
			key = line[0:2].strip()
			data = line[2:].strip()
			result[key] = data

	return result


""" map the data """
model = read_dict('screening_models.txt')
strain = read_dict('strain.txt')
site  = read_dict('site.txt')
schedule = str
route  =   str
vehicle = read_dict('vehicles.txt')
param = read_dict('parameter.txt')
compound = int
unit =  str
dose =  float
tc = int
endpoint = int
survived = str
x = str
mappers = [model, strain, site, schedule, route, vehicle, param, compound, unit, dose, tc, endpoint, survived, x ]



def mapdata(i, value):

	func = mappers[i] 

	if type(func) is dict:
		return func[value] 
	elif func is None:
		return value
	else:
		try:
			return func(value)
		except:
			return value


def duration(schedule):
	#Q01DX005
	try:
		every = int(schedule[1:3])
	except:
		every = 0

	time  = schedule[3]
	try:
		count  = int(schedule[-3:])
	except:
		count = 0

	duration = float(every * count)
	# convert to days
	if time == 'H' : duration /= 24.0
	if time == 'M' : duration /= 24.0 * 60
	return duration


def assess(datamap):
	readout = datamap['parameter']
	result =   datamap['T/C']
	assess = False
	if 'survival time' in readout.lower():
		if result > 130: assess = True
		datamap['active'] = assess
	if 'tumor weight' in readout.lower():
		if result < 20: assess = True
		datamap['active'] = assess
	
	return datamap


def add_tags(sdfile, datamap):
	sdfile = sdfile[:sdfile.find('$$$$')]
	for tag in datamap:
		sdfile += '>  <' + tag + '>\n'
		sdfile += str(datamap[tag]) + '\n\n'

	sdfile += '$$$$'
	return sdfile


def createurl(number):
	queryformat = 'sdf'
	queryurl = ''
	done = False
	while not done:
		try:
			queryurl = 'https://cactus.nci.nih.gov/chemical/structure/NSC%d/file?format=%s' % (number, queryformat)
			done = True
		except:
			time.sleep(5)
			
	return queryurl


def processrecord(i, datamap):

	global cache_hits, cache_miss, cached_structure
	cmpd = int(datamap['NSC'])

	if i > 0 and i % log_interval == 0:
		sys.stderr.write('wrote  %7d sdfiles cache_hit %6d cache_miss %6d \n' % (i, cache_hits, cache_miss) )

	if cmpd in cached_structure:
		cache_hits += 1
		structure = cached_structure[cmpd]
	else:
		cache_miss += 1
		url = createurl(cmpd)
		structure = urlopen(url).read().decode('utf-8')
		cached_structure[cmpd] = structure

	structure = add_tags(structure, datamap)
	print(structure)

#------------------------------start processing --------------------------

# cache structures to lower internet traffics
cached_structure = dict()
cache_hits = 0
cache_miss = 0
log_interval = 1000
lastread = -10

if os.path.exists(ckpoint):
	with open(ckpoint, 'r') as f:
		lastread = int(f.read())
		sys.stderr.write('restarting at line ' +  str(lastread) + '\n')

with open(fname, 'r') as file:
	for count, line in enumerate(file.readlines()):
		if count <= lastread:
			continue

		data = line.split(delim)
		datamap = dict() 
		for i, item in enumerate(data):
			item = item.strip()
			tag  = columns[i]
			item = mapdata(i, item)
			datamap[tag] = item	
	
		if 'schedule' in datamap:	
			datamap['duration'] = duration( datamap['schedule'])

		datamap = assess(datamap)
		processrecord(count, datamap)

		# remember last entry read for restarting
		with open(ckpoint, 'w') as f:
			f.write(str(count))


sys.stderr.write('END')

