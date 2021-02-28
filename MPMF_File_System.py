import os

class FileSystem:
    """
        A structure for the configuration files
        Sets the location of the configuration files
        Used by ProcessRawFile and MPMFDBSetup
        https://realpython.com/documenting-python-code/#documenting-your-python-code-base-using-docstrings
    """
    def __init__(self, file_directory, out_directory, machine, experiment):

        self.main_dir = os.getcwd()
        self.sw_dir = self.main_dir + "\\" + "Software"
        self.config_dir = self.main_dir + "\\" + "Config"
        self.in_dir = file_directory

        # config for processing
        if out_directory != "":
            self.instruments = self.config_dir + "\\" + "instruments"+ "\\" + experiment
            self.xml_template_metab = self.config_dir + "\\" + "metab_template.xml"
            self.xml_template_proteo = self.config_dir + "\\" + "proteo_template.xml"
            self.out_dir = out_directory #"C:\Users\sliggady\Desktop\OutfilesTest"

            # set Proteomics mzmine input file
            if experiment == "proteomics":
                instr_file = "iRT-Reference-" + machine + ".csv"
                self.irt_db = self.instruments + "\\" + instr_file
                if not os.path.exists(self.irt_db):
                    self.irt_db = self.config_dir + "\\" + "iRT-Reference-Default.csv"

                # and threshold file for email notification
                threshold_file = "proteomics-thresholds-" + machine + ".txt"
                self.thresh_email = self.instruments + "\\" + threshold_file
                if not os.path.exists(self.thresh_email):
                    self.thresh_email = self.config_dir + "\\" + "proteomics-thresholds-default.txt"

            # set Metabolomics mzmine input file
            if experiment == "metabolomics":
                instr_file_neg = "negative-db-" + machine + ".csv"
                instr_file_pos = "positive-db-" + machine + ".csv"
                self.neg_db = self.instruments + "\\" + instr_file_neg
                self.pos_db = self.instruments + "\\" + instr_file_pos
                if not os.path.exists(self.neg_db):
                    self.neg_db = self.config_dir + "\\" + "negative-db-Default.csv"
                if not os.path.exists(self.pos_db):
                    self.pos_db = self.config_dir + "\\" + "positive-db-Default.csv"

                # and threshold file for email notification
                threshold_file = "metab-thresholds-" + machine + ".txt"
                self.thresh_email = self.instruments + "\\" + threshold_file
                if not os.path.exists(self.thresh_email):
                    self.thresh_email = self.config_dir + "\\" + "metab-thresholds-default.txt"

if __name__ == "__main__":
    FileSystem("Test", "", "qeclassic", "metabolomics")
