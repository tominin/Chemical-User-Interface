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

def calc_remaining(table, label, all_data):
	current_data = len(table.get_children())   #gets number of rows currently in gui
	text_out = 'Remaining Molecules:\n\n{}/{}'.format(current_data, all_data)    #prints current / when first loaded into gui from sqlite3
	label.configure(text=text_out)   #updates text on label

def tdf_reader(sdf, con=None):
	supply = Chem.SDMolSupplier(sdf)   #rdkit class that reads in sdf and creates an itterable list of moles
	for mol in supply:
		try:
			temp_mw = mol.GetProp('MolWeight')
			mw = temp_mw.replace(',', '')   #if molweight > 1000 includes a ,  in string
			H_acceptors = Descriptors.NumHAcceptors(mol)
			H_donors = Descriptors.NumHDonors(mol)
			num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
			mol_logd = mol.GetProp('LogD')
			temp_logp = Descriptors.MolLogP(mol)
			mol_logp = round(temp_logp, 2)  #round to make prettier
			mol_name = mol.GetProp('_Name')  #gets IUPAC name
			smiles = Chem.MolToSmiles(mol) #creates a smile from mol
			atom_num = mol.GetNumAtoms()
			formula = mol.GetProp("Formula")
			mol_ID = mol.GetProp("Mol_ID")   #used as primark key
			pli = Draw.MolToImage(mol)
			BLOB = pickle.dumps(pli)   #this is the BLOB, for the sql database, just a binary format of the image to be saved and easily converted pack to png when wanted to
			num_rings = rdMolDescriptors.CalcNumRings(mol)
			psa_temp = rdMolDescriptors.CalcTPSA(mol)
			psa = round(psa_temp, 2)
			instance = Molecule_Data(BLOB, mol_ID, formula, mol_name, smiles, atom_num, mw,  mol_logp, mol_logd, H_acceptors, H_donors, num_rotatable_bonds, num_rings, psa)   #sends all the relevant information into Molecule_Date class
			instance.insert_sql(con)
		except AttributeError:
			pass
def sql_reader(output, table):
	for row in output:
		pli = pickle.loads(row[0])   #pickle.loads converts from BLOB to png
		pli = pli.resize((128, 128))
		png = ImageTk.PhotoImage(pli)
		mol_id = row[1].replace(",", "")   #for 1000
		instance = Molecule_Data(png, int(mol_id), row[2], row[3], row[4], int(row[5]), float(row[6]), float(row[7]), float(row[8]),  int(row[9]), int(row[10]), int(row[11]), int(row[12]), float(row[13]))   #sends all relevevant information into Molecule_Date class to be uploaded to gui in same order
		instance.insert_gui(table)   #inserted into gui row by row

def sort_treeview(table, col, descending):   #STACKOVERFLOW/CHATGPT
	data = [(float(table.set(item, col)) if table.set(item, col).replace('.', '', 1).replace('-', '').isdigit() else table.set(item, col), item) for item in table.get_children('')]   #gets the float item for a a selected column
	data.sort(reverse=descending, key=lambda x: (isinstance(x[0], str), x[0]))   #orders by descending value
	for index, (val, item) in enumerate(data):   #this gets a new index position using enumarate to value the list of values in a column by numerical (if float) and this index is used to reposition the row in the gui
		table.move(item, '', index)
	table.heading(col, command=lambda: sort_treeview(table, col, not descending))   #recall function to apply

def text_button_writer(table, widgets, texts):
	for value, widget in zip(texts.values(), widgets):
		widget.set(value)   #relabels the gui entry text boxes after pressing one of the buttons to show the default values
def lipinski_filter(table, removed_table, widgets, texts, all_data, label):
	for item in table.get_children():   #gets the rows of the gui basically
		mw_value = float(table.item(item, "values")[5])
		ha_value = int(table.item(item, "values")[8])
		hd_value = int(table.item(item, "values")[9])
		logP_value = float(table.item(item, "values")[6])
		if mw_value <= 500.0 and ha_value <= 10 and hd_value <= 5 and logP_value <= 5.0:   #list of default parameters
			table.reattach(item, '', 'end')
		else:
			table.detach(item)
			removed_table.append(item)
	texts['mw_value'] = 500   #have to hardcode the values and reset any others back to 0
	texts['ha_value'] = 10
	texts['hd_value'] = 5
	texts['logp_value'] = 5
	texts['logd_value'] = ''   #empty
	texts['psa_value'] = ''  #empty
	texts['rings_value'] = ''
	texts['rot_bonds_value'] = ''
	text_button_writer(table, widgets, texts)  #writes text box
	calc_remaining(table, label, all_data)  #calculates remaining rows
	return (removed_table)  #removed table is the detatched rows to be kept as a global variable to allow for re insertion (reset) button
def lead_likliness_filter(table, removed_table, widgets, texts, all_data, label):    #very similar to lipinski, possibly a better way to do this is in a class and make a superclass to inherit common functions, but the function is the same and I found this easier to code in the time
	for item in table.get_children():
		mw_value = float(table.item(item, "values")[5])
		ha_value = int(table.item(item, "values")[8])
		hd_value = int(table.item(item, "values")[9])
		logD_value = float(table.item(item, "values")[7])
		ring_count_value = int(table.item(item, "values")[11])
		rotatable_bonds_count_value = int(table.item(item, "values")[10])
		if mw_value <= 450 and ha_value <= 8 and hd_value <= 5 and ring_count_value <= 4 and rotatable_bonds_count_value <= 10 and logD_value >= -4 and logD_value <= 4:
			table.reattach(item, '', 'end')
		else:
			table.detach(item)
			removed_table.append(item)
	texts['mw_value'] = 450
	texts['ha_value'] = 8
	texts['hd_value'] = 5
	texts['logp_value'] = ''
	texts['logd_value'] = 4
	texts['psa_value'] = ''  #empty
	texts['rot_bonds_value'] = 10
	texts['rings_value'] = 4
	text_button_writer(table, widgets, texts)
	calc_remaining(table, label, all_data)
	return (removed_table)
def bioavailability_filter(table, removed_table, widgets, texts, all_data, label):
	for item in table.get_children():
		score = 0
		mw_value = float(table.item(item, "values")[5])
		ha_value = int(table.item(item, "values")[8])
		hd_value = int(table.item(item, "values")[9])
		logP_value = float(table.item(item, "values")[6])
		rotatable_bonds_count_value = int(table.item(item, "values")[10])
		ring_count_value = int(table.item(item, "values")[11])
		psa_value = float(table.item(item, "values")[12])
		if mw_value <= 500.0:
			score = score + 1
		if ha_value <= 10:
			score = score + 1
		if hd_value <= 5:
			score = score + 1
		if logP_value <= 5.0:
			score = score + 1
		if rotatable_bonds_count_value <= 10:
			score = score + 1
		if ring_count_value <= 5:
			score = score + 1
		if psa_value <= 200:
			score = score + 1
		if score >= 6:
			table.reattach(item, '', 'end')
		else:
			table.detach(item)
			removed_table.append(item)
	texts['mw_value'] = 500
	texts['ha_value'] = 10
	texts['hd_value'] = 5
	texts['logp_value'] = 5
	texts['logd_value'] = ''  #empty
	texts['rot_bonds_value'] = 10
	texts['rings_value'] = 5
	texts['psa_value'] = 200
	text_button_writer(table, widgets, texts)
	calc_remaining(table, label, all_data)
	return (removed_table)

def reset_table(table, removed_table, all_data, label):
	for item in removed_table:
		table.reattach(item, '', 'end')   #reattatch tk function is amazing this is such a nice easy bit of code for what it does
	calc_remaining(table, label, all_data)

def apply_new_filters(table, removed_table, all_data, label, widgets):
	fixed_criteria = [float(x.get()) if x.get().replace('.', '').replace('-', '').isdigit() else 9999999 for x in widgets]    #sets values as float but if a misentry or nothing set to 9999999 (as all filters filter as a max) this would be something to improve on a better version possibly having a min and max
	for item in table.get_children():
		mw_value = float(table.item(item, "values")[5])
		logP_value = float(table.item(item, "values")[6])
		logD_value = float(table.item(item, "values")[7])
		psa_value = float(table.item(item, "values")[12])
		ha_value = int(table.item(item, "values")[8])
		hd_value = int(table.item(item, "values")[9])
		ring_count_value = int(table.item(item, "values")[11])
		rotatable_bonds_count_value = int(table.item(item, "values")[10])
		if mw_value <= fixed_criteria[0] and logP_value <= fixed_criteria[1] and logD_value <= fixed_criteria[2] and logD_value >= -(fixed_criteria[2]) and psa_value <= fixed_criteria[3] and ha_value <= fixed_criteria[4] and hd_value <= fixed_criteria[5] and ring_count_value <= fixed_criteria[6] and rotatable_bonds_count_value <= fixed_criteria[7]:
			table.reattach(item, '', 'end')
		else:
			table.detach(item)
			removed_table.append(item)
	calc_remaining(table, label, all_data)
	return (removed_table)

def extract_table(table):   #would of liked to include more information in SDF file in a more refined programme.
	with Chem.SDWriter("shortlisted.sdf") as writer:
		for item in table.get_children():
			smiles = str(table.item(item, "values")[3])
			if smiles != None:
				mol = Chem.MolFromSmiles(smiles)
				writer.write(mol)

class Molecule_Data:
	def __init__(self, BLOB, mol_ID, formula, IUPAC, SMILES, atom_number, mw, logP, logD, HA, HD, num_rotatable_bonds, num_rings, psa):   #class of question
		self._mol_ID = mol_ID
		self._psa = psa
		self._num_rings = num_rings
		self._BLOB = BLOB
		self._formula = formula
		self._IUPAC = IUPAC
		self._SMILES = SMILES
		self._atom_number = int(atom_number)
		self._mw = float(mw)
		self._logP = float(logP)
		self._logD = float(logD)
		self._HA = int(HA)
		self._HD = int(HD)
		self._num_rotatable_bonds = int(num_rotatable_bonds)
	@property
	def psa(self):
		return self._psa
	@property
	def num_rings(self):
		return self._num_rings
	@property
	def logD(self):
		return self._logD
	@property
	def mol_ID(self):
		return self._mol_ID
	@property
	def BLOB(self):
		return self._BLOB
	@property
	def formula(self):
		return self._formula
	@property
	def IUPAC(self):
		return self._IUPAC
	@property
	def SMILES(self):
		return self._SMILES
	@property
	def atom_number(self):
		return self._atom_number
	@property
	def mw(self):
		return self._mw
	@property
	def logP(self):
		return self._logP
	@property
	def HA(self):
		return self._HA
	@property
	def HD(self):
		return self._HD
	@property
	def num_rotatable_bonds(self):
		return self._num_rotatable_bonds
	def insert_sql(self, con):
		sql_i = "INSERT INTO Molecule_Data VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
		params_i = (self.BLOB, self.mol_ID, self.formula, self.IUPAC, self.SMILES, self.atom_number, self.mw, self.logP, self.logD, self.HA, self.HD, self.num_rotatable_bonds, self.num_rings, self.psa)
		try:
			cur = con.cursor()
			rows = cur.execute(sql_i, params_i)
			con.commit()
			cur.close()
		except:
			print ("Entry {} already inserted".format(self.mol_ID))
	def insert_gui(self, table):
		if not hasattr(table, 'image_refs'):
			table.image_refs= []
		table.insert(parent="", index="end", text="", image=self.BLOB, values=(self.mol_ID, self.formula, self.IUPAC, self.SMILES, self.atom_number, self.mw, self.logP, self.logD, self.HA, self.HD, self.num_rotatable_bonds, self.num_rings, self.psa), tags=('rowstyle', ))
		table.image_refs.append(self.BLOB)
