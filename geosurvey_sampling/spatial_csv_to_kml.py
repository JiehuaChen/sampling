#!/usr/bin/python

import os, sys, glob, shutil, getopt, re
import csv

def main():				
	csv_files = glob.glob("csv/*.csv")
	for f in csv_files:
		print "Writing KML file for %s" % (f,)
		csv_to_kml(f)	
	print "Done."

def csv_to_kml(input_filename):

	# open input file
	csv_file = open(input_filename,'rU')
	reader = csv.DictReader(csv_file)
	# preamble 
	input_filename_base, input_filename_ext = os.path.splitext(input_filename)
		
	# open output file
	kml_file = open(input_filename_base + '.kml','w')

	kml_file.write(r"""<?xml version="1.0" encoding="utf-8" ?>
	<kml xmlns="http://www.opengis.net/kml/2.2">
	""")
	kml_file.write("<Document><name>%s</name>" % input_filename_base)
	kml_file.write(r"""
	<Schema name="sample" id="sample">
		<SimpleField name="Name" type="string"></SimpleField>
		<SimpleField name="Description" type="string"></SimpleField>
		<SimpleField name="L2" type="int"></SimpleField>
		<SimpleField name="L1" type="int"></SimpleField>
		<SimpleField name="L0" type="int"></SimpleField>
		<SimpleField name="ID_unit" type="int"></SimpleField>
		<SimpleField name="Prob_3_stage" type="float"></SimpleField>
		<SimpleField name="Prob" type="float"></SimpleField>
	</Schema>
	""")

	# main loop 
	for line in reader:

		kml_file.write('  <Placemark>\n')
		kml_file.write('  <name>L2=%s, L1=%s, L0=%s</name>\n' % (line['L2'],line['L1'],line['L0']))
		kml_file.write('\t<ExtendedData><SchemaData schemaUrl=\"#sample\">\n')
		kml_file.write(' <SimpleField name="ID_unit">%s</SimpleField>\n' % (line['ID_unit']))
		kml_file.write('  <SimpleField name="Prob_3_stage">%s</SimpleField>\n' % (line['Prob_ 3 _stage']))
		kml_file.write('  <SimpleField name="Prob">%s</SimpleField>\n' % (line['Prob']))
		kml_file.write('\t\t</SchemaData></ExtendedData>\n')
		kml_file.write("     <Point><coordinates>%s,%s</coordinates></Point>\n" % (line['x'], line['y']))
		kml_file.write('  </Placemark>\n')
	
	# epilogue
	kml_file.write('\t</Document>\n\t</kml>')
	csv_file.close()
	kml_file.close()

if __name__ == '__main__':
	main()
