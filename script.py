#!/home/tom/Documents/glasgow_uni/semester2/chemical_databases/python_db_assessment/chem_db_env/bin/python3
import argparse
import sqlite3
import rdkit
import re
import logging
import os
import webbrowser
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
import pickle
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
import io
import ChemDb as CDb

##
parser = argparse.ArgumentParser(description='args for file inputs', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--sdf_file', help='Required entry, enter the full file path for the sdf file you want to check')
parser.add_argument('--database', help='Required entry, enter the full file path for sql database file')
args = parser.parse_args()
sdf_file = args.sdf_file
db_file = args.database
##

if not os.path.isfile(db_file):   #if database argument is created then it just connects and carries on but if it is not created will create and connect
	print ('Creating database...\n')
	with open(db_file, 'w') as o:
		con = sqlite3.connect(db_file)
	print ('Successfully created empty database in current directory called "omics.db"\n')
	try:
##creates a sqlite3 schema making one table of molecule data for everything later required
		con.executescript("""
BEGIN;
CREATE TABLE Molecule_Data(
		Molecular_Structure BLOB,
		MolID varchar(9) NOT NULL,
		Formula varchar(50),
		IUPAC varchar(25),
		SMILES varchar(25) ,
		Atom_Number int,
		Molecular_Weight decimal(12, 6),
		logP decimal(8, 6),
		logD decimal(8, 6),
		HA int,
		HD int,
		Number_of_Rotatable_Bonds int,
		Ring_Count int,
		PSA decimal(12, 6),
		PRIMARY KEY(MolID)
);
COMMIT;
""")
	except sqlite3.OperationalError:
		print ("database already made")   #this in reality would never be triggered, but just here in case!
	CDb.tdf_reader(sdf_file, con)   #sends the sdf_file to a class that reads in all the data and uploads it to the database line by line, changing them into correct data types

else:
	print ('Connecting to existing database...\n')
	con = sqlite3.connect(db_file)  #if database is already connected then it will skip to straight up taking all the data out of the sqlite3 database

app = tk.Tk()   ##app = root window
app.title("BioActive Molecules Database")   #name of the gui
columns = ("MolID", "Formula", "IUPAC", "SMILES", "Atom Number", "Molecular Weight", "logP", "logD", "Hydrogen Acceptors", "Hydrogen Donors", "Number of Rotatable Bonds", "Number of Rings", "PSA")  #all cols of gui
table = ttk.Treeview(app, columns=columns)   #table is the master root
for col in columns:
	table.heading(col, text=col, command=lambda c=col: CDb.sort_treeview(table, c, False))   #this calls sort_treeview a function that just orders the gui treeview by value when selected

table.column('#0', width=140)   #first row holds the images
table.column('MolID', anchor=tk.W, width=20)
table.column('Formula', anchor=tk.W, width=100)
table.column('IUPAC', anchor=tk.W, width=100)
table.column('SMILES', anchor=tk.W, width=200)
table.column('Atom Number', anchor=tk.W, width=50)
table.column('Molecular Weight', anchor=tk.W, width=80)
table.column('logP', anchor=tk.W, width=50)
table.column('logD', anchor=tk.W, width=50)
table.column('Hydrogen Acceptors', anchor=tk.W, width=50)
table.column('Hydrogen Donors', anchor=tk.W, width=50)
table.column('Number of Rotatable Bonds', anchor=tk.W, width=50)
table.column('Number of Rings', anchor=tk.W, width=50)
table.column('PSA', anchor=tk.W, width=50)
table.heading("#0", text='Molecular Structure')   #labels of the columns below
table.heading('MolID', text='MolID', anchor=tk.W)
table.heading('Formula', text='Formula', anchor=tk.W)
table.heading('IUPAC', text='IUPAC', anchor=tk.W)
table.heading('SMILES', text='SMILES', anchor=tk.W)
table.heading('Atom Number', text='Atom Number', anchor=tk.W)
table.heading('Molecular Weight', text='Molecular Weight', anchor=tk.W)
table.heading('logP', text='logP', anchor=tk.W)
table.heading('logD', text='logD', anchor=tk.W)
table.heading('Hydrogen Acceptors', text='Hydrogen Acceptors', anchor=tk.W)
table.heading('Hydrogen Donors', text='Hydrogen Donors', anchor=tk.W)
table.heading('Number of Rotatable Bonds', text='Number of Rotatable Bonds', anchor=tk.W)
table.heading('Number of Rings', text='Number of Rings', anchor=tk.W)
table.heading('PSA', text='Topological PSA', anchor=tk.W)
table.tag_configure('rowstyle', background='white smoke')   #all rows coloured in white smoke

cur = con.cursor()
cur.execute('SELECT * FROM Molecule_Data')   #takes all the data from sqlite3 database
output = cur.fetchall()
all_data = len(output)   #stored for later use
CDb.sql_reader(output, table)   #class to read in the data from this and then upload to gui
cur.close()

removed_table = []
##LEFT BUTTON COLUMN
button_frame = tk.Frame(app, bg="thistle3")   #creating a seperate frame for the buttons on the left side
button_frame.pack(side="left", fill=tk.BOTH)
criteria_label = tk.Label(button_frame, text="Personalised Criteria", bg="thistle3")   #just to make the gui look clearer
criteria_label.pack(padx=5, pady = (25, 50))
lipinski_button = tk.Button(button_frame, text="Lipinksi's Filter", bg="purple1", command= lambda: CDb.lipinski_filter(table, removed_table, widgets, texts, all_data, current_label)) #when pressed filters for lipinski
lead_drug_button = tk.Button(button_frame, text="Lead-likeness Filter", bg="purple1", command= lambda: CDb.lead_likliness_filter(table, removed_table, widgets, texts, all_data, current_label))  #filters for lead-drug
bio_button = tk.Button(button_frame, text="Bioavailability Filter", bg="purple1", command= lambda: CDb.bioavailability_filter(table, removed_table, widgets, texts, all_data, current_label))  #filters for bioavailability
reset_button = tk.Button(button_frame, text="Reset", bg="purple1", command= lambda: CDb.reset_table(table, removed_table, all_data, current_label)) #restores dataase to beginning
extract_button = tk.Button(button_frame, text="Extract to SDF", bg="purple1", command = lambda: CDb.extract_table(table))  #button that sends all current entries into a sdf file
current_label = tk.Label(button_frame, text="", bg="thistle3")   #this is initially called as empty but a function taht is called later and everytime any of the buttons are pressed will tell the user how many of the original molecules are still present after filters applied
lipinski_button.pack(padx=5, pady=(100,50))
lead_drug_button.pack(padx=5, pady=50)
bio_button.pack(padx=5, pady=50)
reset_button.pack(padx=5, pady=50)
extract_button.pack(padx=5, pady=50)
current_label.pack(padx=5, pady=50)
CDb.calc_remaining(table, current_label, all_data)  #this the function in question

##TOP TEXT INPUT ROWS
##first button row
text_frame = tk.Frame(app, bg="thistle3")
text_frame.pack(side="top", fill=tk.BOTH)
mw_text_label = tk.Label(text_frame, text="MW:", bg="thistle3")
mw_text_label.grid(row = 0, column = 0, padx = 50)
mw_entry = tk.StringVar()
mw_field = tk.Entry(text_frame, textvariable=mw_entry)
mw_field.grid(row = 0, column = 1, padx = 50)
logp_text_label = tk.Label(text_frame, text="LogP:", bg="thistle3")
logp_text_label.grid(row = 0, column = 2, padx = 50)
logp_entry = tk.StringVar()
logp_field = tk.Entry(text_frame, textvariable=logp_entry)
logp_field.grid(row = 0, column = 3, padx = 50)
logd_text_label = tk.Label(text_frame, text="LogD:", bg="thistle3")
logd_text_label.grid(row = 0, column = 4, padx = 50)
logd_entry = tk.StringVar()
logd_field = tk.Entry(text_frame, textvariable=logd_entry)
logd_field.grid(row = 0, column = 5, padx = 50)
psa_text_label = tk.Label(text_frame, text="PSA:", bg="thistle3")
psa_text_label.grid(row = 0, column = 6, padx = 50)
psa_entry = tk.StringVar()
psa_field = tk.Entry(text_frame, textvariable=psa_entry)
psa_field.grid(row = 0, column =7, padx = 50)
##second button row
ha_text_label = tk.Label(text_frame, text="HA:", bg="thistle3")
ha_text_label.grid(row = 1, column = 0, padx = 50, pady = 10)
ha_entry = tk.StringVar()
ha_field = tk.Entry(text_frame, textvariable=ha_entry)
ha_field.grid(row = 1, column =1, padx = 50, pady = 10)
hd_text_label = tk.Label(text_frame, text="HD:", bg="thistle3")
hd_text_label.grid(row = 1, column = 2, padx = 50, pady = 10)
hd_entry = tk.StringVar()
hd_field = tk.Entry(text_frame, textvariable = hd_entry)
hd_field.grid(row = 1, column =3, padx = 50, pady = 10)
rings_text_label = tk.Label(text_frame, text="Num of Rings:", bg="thistle3")
rings_text_label.grid(row = 1, column = 4, padx = 50, pady = 10)
rings_entry = tk.StringVar()
rings_field = tk.Entry(text_frame, textvariable = rings_entry)
rings_field.grid(row = 1, column =5, padx = 50, pady = 10)
rot_bonds_text_label = tk.Label(text_frame, text="Num of Rotatable Bonds:", bg="thistle3")
rot_bonds_text_label.grid(row = 1, column = 6, padx = 10, pady = 10)
rot_bonds_entry = tk.StringVar()
rot_bonds_field = tk.Entry(text_frame, textvariable = rot_bonds_entry)
rot_bonds_field.grid(row = 1, column =7, padx = 10, pady = 10)

texts = {'mw_value': None, 'logp_value': None, 'logd_value': None, 'psa_value': None, 'ha_value': None, 'hd_value': None, 'rings_value': None, 'rot_bonds_value': None}   #dict of the customisable variables used to be sent into all the functions to update
widgets = [mw_entry,logp_entry,logd_entry,psa_entry,ha_entry,hd_entry,rings_entry,rot_bonds_entry]   #list of all entries, important taht order is the same as the dict
submit_button = tk.Button(text_frame, text="Submit", command = lambda: CDb.apply_new_filters(table, removed_table, all_data, current_label, widgets), bg="purple1")   #submit press to apply customisable filters
submit_button.grid(row = 3, column = 4, pady = (0, 10))

def print_select(e):
	sel = table.focus()   #.focus sends the highlighted row data to sel variable
	try:
		name = str(table.item(sel, "values")[2])
		webbrowser.open("https://pubchem.ncbi.nlm.nih.gov/compound/{}".format(name), new=1)   #opens pubchem

	except:
		pass
def big_image(e):
	sel = table.focus()
	try:
		mol_id = table.item(sel, "values")[0]
		cur = con.cursor()
		cur.execute('SELECT * FROM Molecule_Data WHERE MolID IS {}'.format(mol_id))
		output = cur.fetchone()
		cur.close()
		pli = pickle.loads(output[0])
		pli = pli.resize((600, 600))
		png = ImageTk.PhotoImage(pli)
		novi = tk.Toplevel(table)
		novi.geometry("600x600")   #biggere image
		label = tk.Label(novi, image=png)   #new window!
		label.image = png
		label.pack()
	except:
		pass

table.pack(expand=True, padx=(0, 0), fill=tk.BOTH)
table.bind("<Double-1>", print_select)   #double left click to see the molecules pubchem page
table.bind("<Double-3>", big_image)   #double right click to get an enlarged image
vsb = ttk.Scrollbar(app, orient="vertical", command=table.yview)   #vertical scroll  bar
vsb.pack(side="right", fill="y")
table.configure(yscrollcommand=vsb.set)
hsb=ttk.Scrollbar(app, orient="horizontal", command = table.xview)   #horizontal scroll bar
hsb.pack(side="bottom", fill="x")
table.configure(xscrollcommand=hsb.set)
style = ttk.Style()
style.configure("Treeview", rowheight=140)   #rowheight is set this size to nicely fit in the mol images
app.mainloop()   #opens the gui
