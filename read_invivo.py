from urllib.request import urlopen
import pubchempy
from pubchempy import get_compounds
from rdkit import Chem
import sys
import os
import time

fname = 'invivodec02-sorted' # file name
output_file = 'NCI-data.sdf'

ckpoint = 'checkpoint.cpk'

USE_PUBMED = False

delim = '|'  # file delimiter

# columns in file.  I don't know what the last one is
columns = ['tumor', 'strain', 'site', 'schedule', 'route', 'vehicle', 'parameter', 'NSC', 'unit', 'dose', 'T/C', 'endpoint', 'survived', 'x' ]

def errormsg(msg):
	sys.stderr.write(msg + '\n')


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
	if time == 'M' : duration /= (24.0 * 60)
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
	if sdfile[-1] != '\n' : sdfile += '\n'

	for tag in datamap:
		sdfile += '>  <' + tag + '>\n'
		sdfile += str(datamap[tag]) + '\n\n'

	sdfile += '$$$$\n'
	return sdfile


def getsdf(identifier):
	if USE_PUBMED:
		return getsdf_pubmed(identifier)
	else:
		return getsdf_nci(identifier)


def getsdf_nci(number):
	queryformat = 'sdf'
	timeout = 10
	sdfile = None
	number = int(number)
	queryurl = 'https://cactus.nci.nih.gov/chemical/structure/NSC%d/file?format=%s' % (number, queryformat)

	for i in range(2):
		try:
			sdfile = urlopen(queryurl,timeout=timeout).read().decode('utf-8')

		except Exception as e:
			errormsg('sleeping: ' + str(e))
			time.sleep(5)

	if sdfile == None:	
		errormsg('error fetching ' + queryurl  + ' ' + str(e))
		errormsg('trying pubmed')
		sdfile = getsdf_pubmed(number)

	return sdfile


def getsdf_pubmed(identifier):
	sdfil = None

	# try a few times in case is busy
	for i in range(5):
		try:
			compound = get_compounds('NSC' + str(identifier), 'name')
			break
		except:
			errormsg('sleeping')
			time.sleep(10)


	if compound is not None and len(compound) > 0:
		compound = compound[0]
		try:
			smiles = compound.isomeric_smiles 
			mol = Chem.MolFromSmiles(smiles)
			sdfil = Chem.MolToMolBlock(mol)
		except Exception as e:
			errormsg('error converting smiles to sdf: ' + smiles)
			errormsg(str(e))
			

	return sdfil	

def processrecord(datamap):

	global cache_hits, cache_miss, cached_structure
	structure = None
	cmpd = int(datamap['NSC'])


	if cmpd in cached_structure:
		cache_hits += 1
		structure = cached_structure[cmpd]
	else:
		cache_miss += 1
		structure = getsdf(cmpd)
		if len(cached_structure) > max_cache:
			cached_structure = dict()

		cached_structure[cmpd] = structure

	if structure is not None:
		structure = add_tags(structure, datamap)
		with open(output_file, 'a') as f:
			f.write(structure) # has newline at end already
		return 0
	else:
		errormsg('error finding ' + str(cmpd) )
		return 1

#------------------------------start processing --------------------------

# cache structures to lower internet traffics
cached_structure = dict()
cache_hits = 0
cache_miss = 0
max_cache = 100
log_interval = 100
lastread = -10

if os.path.exists(ckpoint):
	with open(ckpoint, 'r') as f:
		lastread = int(f.read())
		errormsg('restarting at line ' +  str(lastread))

with open(fname, 'r') as file:
	start = time.time()
	errorcount = 0
	for count, line in enumerate(file.readlines()):
		if count <= lastread:
			continue

		if count  > 0 and count  % log_interval == 0:
			elapsed = time.time() - start
			recs_sec = count/elapsed
			errormsg('wrote  %7d sdfiles cache_hit %6d cache_miss %6d error %5d rate %5.2f' 
				% (count, cache_hits, cache_miss, errorcount, recs_sec) )

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
		errorcount += processrecord(datamap)

		# remember last entry read for restarting
		with open(ckpoint, 'w') as f:
			f.write(str(count))

errormsg('END')

