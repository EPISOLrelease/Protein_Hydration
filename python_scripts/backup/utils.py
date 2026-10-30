import numpy as np
import os
class Molecule(dict):
    def __init__(self, file):
        from _io import TextIOWrapper
        read_type = False

        if (type(file) == TextIOWrapper):
            self.filetype = file.name[-3:]
        else:
            self.filetype = file[-3:]
            read_type = True
        #probably a bad way of doing that
        #print(self.filetype)
        file_type_ = self.filetype
        ## KEYS
        self["REMARK"] = ""
        self[f"{file_type_}"] = []
        ##########################
        # okay need to make this a @classmethod
        if read_type:
            with open(file,'r') as rr:
                if file_type_ == 'pdb':
                    self.pdb_reader(rr)
                if file_type_ == 'gro':
                    self.gro_reader(rr)
        #probably a bad way of doing that
        else:
            if file_type_ == 'pdb':
                self.pdb_reader(file)
            if file_type_ == 'gro':
                self.gro_reader(file)
        ##########################
    def __str__(self):
        outstr = ""
        for at in self:
            outstr += str(at)

        return outstr
    ##############################
    # FUNCTIONS
    ##############################
    #@classmethod
    def pdb_reader(cls,file):
        file_type_ = cls.filetype
        for line_count ,line in enumerate(file):
            try:
                if line[:6] != "REMARK":
                    if "ATOM" in line or "HETATM" in line:
                        # list ver self.append(Atom(line,self.filetype))
                        cls[f"{file_type_}"].append(Atom(line,file_type_))
            except KeyboardInterrupt:
                return  
            except Exception as exc:
                #print(f"ERROR ON LINE {line_count} OF INPUT FILE {self.file}")
                print(f"ERROR: {RuntimeWarning(exc)} AT INPUT FILE LINE: {line_count}")
                raise RuntimeWarning(exc)
    #@classmethod
    def gro_reader(cls,file):
        file_type_ = cls.filetype
        atom_count = int()
        for line_count ,line in enumerate(file):
            #if line_count < 1: continue
            try:
                if line_count == 1: cls.num_atoms = atom_count = int(line.split()[0])
                if line_count == atom_count+2: break
                if line_count > 1: # ignore header
                    #if "ATOM" in line or "HETATM" in line:
                    # list ver self.append(Atom(line,self.filetype))
                    cls[f"{file_type_}"].append(Atom(line,file_type_))
            except KeyboardInterrupt:
                return  
            except Exception as exc:
                #print(f"ERROR ON LINE {line_count} OF INPUT FILE {self.file}")
                print(f"ERROR: {RuntimeWarning(exc)} AT INPUT FILE LINE: {line_count}")
                raise RuntimeWarning(exc)
## NOW WE HAVE THE ATOM CLASS (ALSO A DICT)
# doesnt need to inherit anything
class Atom(dict):
    #from numpy import array
    def __init__(self,line,filetype):
        #from numpy import array
        if filetype == 'pdb':
            self["type"] = line[0:6].strip()
            self["idx"] = line[6:11].strip()
            self["name"] = line[12:16].strip()
            self["resname"] = line[17:20].strip()
            self["resid"] = int(int(line[22:26]))
            self["x"] = float(line[30:38])
            self["y"] = float(line[38:46])
            self["z"] = float(line[46:54])
            self["occupancy"] = float(line[54:60])
            self["b"] = float(line[60:66])
            self["sym"] = line[76:78].strip()
            ## set params 
            self.type = self["type"]
            self.idx = self["idx"]
            self.name = self["name"]
            self.resname = self["resname"]
            self.resid = self["resid"]
            self.coords = np.array([self["x"],self["y"],self["z"]])
            self.b = self["b"]
            self.occupancy = self["occupancy"]
            self.sym = self["sym"]
        #### END PDB
        elif filetype == 'gro':
            self["resid"] = line[:5].strip()#[0]
            self["resname"] = line[5:10].strip()#[0]
            self["name"] = line[10:15].strip()
            self["idx"] = line[15:20].strip()
            self["x"] = float(line[20:28])
            self["y"] = float(line[28:36])
            self["z"] = float(line[36:44])
            
            ## set params 
            #self.type = self["type"]
            self.idx = self["idx"]
            self.name = self["name"]
            self.resname = self["resname"]
            self.resid = self["resid"]
            self.coords = np.array([self["x"],self["y"],self["z"]])
            
    ## PRINT ONLY WORKS WITH PDB    
    def __str__(self):
        line = list(" " * 80)

        line[0:6] = self["type"].ljust(6)
        line[6:11] = self["idx"].ljust(5)
        line[12:16] = self["name"].ljust(4)
        line[17:20] = self["resname"].ljust(3)
        line[22:26] = str(self["resid"]).ljust(4)
        line[30:38] = f"{self['x']:8.3f}".rjust(8)
        line[38:46] = f"{self['y']:8.3f}".rjust(8)
        line[46:54] = f"{self['z']:8.3f}".rjust(8)
        line[54:60] = f"{self['occupancy']:6.2f}".rjust(6)
        line[60:66] = str(self['b']).rjust(6)
        line[76:78] = self['sym'].rjust(2)
        return "".join(line) + "\n"
