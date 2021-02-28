from MPMF_File_System import FileSystem

try:
    import pymysql as MySQLdb
    MySQLdb.install_as_MySQLdb()
except ImportError:
    pass
import logging


class MPMFDBSetUp:
    """
        Database access module
        When run as a script sets up a DB based on config files
        Create a database in MYSQL before running
        Refer to QC ER Diagram for design
        https://realpython.com/documenting-python-code/#documenting-your-python-code-base-using-docstrings
    """

    # CONSTRUCTOR connects to database as user with pword
    def __init__(self, user, pword, database, filesystem):
        self.username = user
        self.password = pword
        self.database = database
        self.fs = filesystem
        self.create_logger()
        self.connected = False
        try:
            self.db = MySQLdb.connect(host="localhost", user=self.username, password=self.password, db=self.database)
            self.cursor = self.db.cursor()
            #self.logger.warning("Database Connection Success")
            self.connected = True
        except Exception as e:
            self.logger.exception(e)

    def set_up(self):
        self.drop_all_tables()
        self.create_all_tables()
        self.insert_all()
        self.select_and_display_all('metric')
        self.select_and_display_all('sample_component')
        self.select_and_display_all('machine')
        self.select_and_display_all('experiment')

    def create_logger(self):

        # create class level logger as used in processing
        self.logger = logging.getLogger("DATABASE")

        # add handlers for file and console
        c_handler = logging.StreamHandler()
        f_handler = logging.FileHandler('database.log')
        c_handler.setLevel(logging.DEBUG)
        f_handler.setLevel(logging.INFO)

        # add formatting
        c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        c_handler.setFormatter(c_format)
        f_handler.setFormatter(f_format)

        # Add handlers to the logger
        self.logger.addHandler(c_handler)
        self.logger.addHandler(f_handler)
    
    # CREATE TABLES
    def create_all_tables(self):
        self.create_table_experiment()
        self.create_table_machine()
        self.create_table_sample_component()
        self.create_table_metric()
        self.create_table_qc_run()
        self.create_table_measurement()
        self.create_table_stat()

    def create_table_sample_component(self):
        sql = "CREATE TABLE IF NOT EXISTS sample_component (" \
		"component_id INT AUTO_INCREMENT NOT NULL," \
		"component_name TINYTEXT NOT NULL," \
		"component_description TEXT," \
        "component_mode CHAR(1)," \
        "exp_mass_charge DECIMAL(8, 5) NOT NULL," \
        "exp_rt DECIMAL(5, 2) NOT NULL," \
        "experiment_id INT NOT NULL, " \
        "PRIMARY KEY(component_id)," \
        "FOREIGN KEY (experiment_id) REFERENCES experiment(experiment_id)) "
                
        try:
            self.cursor.execute(sql)
        except Exception as e:
            self.logger.exception(e)
            
    def create_table_qc_run(self):
        # no id to avoid accidental duplicates ??
        sql = "CREATE TABLE IF NOT EXISTS qc_run(" \
		"run_id INT AUTO_INCREMENT NOT NULL," \
		"file_name TINYTEXT NOT NULL," \
		"date_time DATETIME NOT NULL," \
        "machine_id INT NOT NULL," \
        "experiment_id INT NOT NULL,"\
        "completed VARCHAR(1) NOT NULL," \
        "PRIMARY KEY(run_id)," \
        "FOREIGN KEY (machine_id) REFERENCES machine(machine_id)," \
        "FOREIGN KEY (experiment_id) REFERENCES experiment(experiment_id))"

        try:
            self.cursor.execute(sql)
        except Exception as e:
            self.logger.exception(e)
            
    def create_table_metric(self):
        sql = "CREATE TABLE IF NOT EXISTS metric (" \
		"metric_id INT AUTO_INCREMENT NOT NULL," \
		"metric_name TINYTEXT," \
		"metric_description TEXT,"  \
        "display_order TINYINT," \
        "display_name TINYTEXT," \
        "use_metab VARCHAR(1) NOT NULL," \
        "use_prot VARCHAR(1) NOT NULL," \
		"morpheus VARCHAR(1) NOT NULL," \
        "PRIMARY KEY(metric_id))"
                
        try:
            self.cursor.execute(sql)
        except Exception as e:
            self.logger.exception(e)
            
    def create_table_measurement(self):
        # watch FLOAT for value for range, powers in area metric
        
        sql = "CREATE TABLE IF NOT EXISTS measurement (" \
		"metric_id INT NOT NULL," \
		"component_id INT NOT NULL," \
		"run_id INT NOT NULL," \
		"value DECIMAL(36, 18)," \
		"PRIMARY KEY(metric_id, component_id, run_id)," \
		"FOREIGN KEY (metric_id) REFERENCES metric(metric_id)," \
		"FOREIGN KEY (component_id) REFERENCES sample_component(component_id)," \
		"FOREIGN KEY (run_id) REFERENCES qc_run(run_id) ON DELETE CASCADE)"
                
        try:
            self.cursor.execute(sql)
        except Exception as e:
            self.logger.exception(e)

    def create_table_stat(self):

        sql = "CREATE TABLE IF NOT EXISTS stat (" \
              "metric_id INT NOT NULL," \
              "component_id INT NOT NULL," \
              "machine_id INT NOT NULL," \
              "count INT," \
              "mean DECIMAL(36, 18)," \
              "std DECIMAL(36, 18)," \
              "min DECIMAL(36, 18)," \
              "25_percent DECIMAL(36, 18)," \
              "50_percent DECIMAL(36, 18)," \
              "75_percent DECIMAL(36, 18)," \
              "max DECIMAL(36, 18)," \
              "PRIMARY KEY(metric_id, component_id, machine_id)," \
              "FOREIGN KEY (metric_id) REFERENCES metric(metric_id)," \
              "FOREIGN KEY (component_id) REFERENCES sample_component(component_id)," \
              "FOREIGN KEY (machine_id) REFERENCES machine(machine_id))"

        try:
            self.cursor.execute(sql)
        except Exception as e:
            self.logger.exception(e)

    def create_table_machine(self):

        sql = "CREATE TABLE IF NOT EXISTS machine (" \
              "machine_id INT AUTO_INCREMENT NOT NULL," \
              "machine_name TINYTEXT NOT NULL," \
              "machine_serial TEXT," \
              "machine_description TEXT," \
              "machine_venue TEXT NOT NULL," \
              "use_metab VARCHAR(1) NOT NULL," \
              "use_prot VARCHAR(1) NOT NULL," \
              "PRIMARY KEY(machine_id))"

        try:
            self.cursor.execute(sql)
        except Exception as e:
            self.logger.exception(e)

    def create_table_experiment(self):
        sql = "CREATE TABLE IF NOT EXISTS experiment (" \
              "experiment_id INT AUTO_INCREMENT NOT NULL," \
              "experiment_type TEXT NOT NULL," \
              "PRIMARY KEY(experiment_id))"

        try:
            self.cursor.execute(sql)
        except Exception as e:
            self.logger.exception(e)
    
    # DROP TABLES
    def drop_table(self, tablename):
        sql = "DROP TABLE " + tablename
        
        try:
            self.cursor.execute(sql)
        except Exception as e:
            self.logger.exception(e)
        self.db.commit()
            
    def drop_all_tables(self):
        tables = ['measurement', 'stat','qc_run', 'sample_component', 'metric', 'machine', 'experiment']
        for table in tables:
            self.drop_table(table)

    # INSERTS
    def insert_all(self):
        self.insert_experiments()
        self.insert_components() 
        self.insert_mzmine_metrics()
        self.insert_morpheus_metrics()
        self.insert_machines()

    def insert_components(self):
        # NOTE: exp_rt's are the default values
        #       can be removed from db

        # get experiment type ids
        sql = "SELECT experiment_id FROM experiment where experiment_type = 'metabolomics'"
        try:
            self.cursor.execute(sql)
            metab_exp_id = self.cursor.fetchone()[0]
        except Exception as e:
            self.logger.exception(e)

        sql = "SELECT experiment_id FROM experiment where experiment_type = 'proteomics'"
        try:
            self.cursor.execute(sql)
            prot_exp_id = self.cursor.fetchone()[0]
        except Exception as e:
            self.logger.exception(e)

        # insert neg metab components
        with open(self.fs.config_dir + "\\" + "negative-db-Default.csv", 'r') as infile:
            for line in infile:
                component = line.strip().split("|")
                if component[0] != 'mz':
                    sql = "INSERT INTO sample_component(component_id, component_name, component_mode, exp_mass_charge, exp_rt, experiment_id) " \
                          "VALUES (NULL, " + "'" + component[2].strip() + "'," + "'N'," + "'" + component[0].strip() + "'," + "'" + component[1].strip() + "','" \
                          + str(metab_exp_id) + "')"

                    try:
                        self.cursor.execute(sql)
                    except Exception as e:
                        self.logger.exception(e)

            self.db.commit()

        # insert pos metab comoponents
        with open(self.fs.config_dir + "\\" + "positive-db-Default.csv", 'r') as infile:
            for line in infile:
                component = line.strip().split("|")
                if component[0] != 'mz':
                    sql = "INSERT INTO sample_component(component_id, component_name, component_mode, exp_mass_charge, exp_rt, experiment_id) " \
                          "VALUES (NULL, " + "'" + component[2].strip() + "'," + "'P'," + "'" + component[0].strip() + "'," + "'" + component[1].strip() + "','" \
                          + str(metab_exp_id) + "')"

                    try:
                        self.cursor.execute(sql)
                    except Exception as e:
                        self.logger.exception(e)

            self.db.commit()

        # insert iRT peptides
        with open(self.fs.config_dir + "\\" + "iRT-Reference-Default.csv", 'r') as infile:
            for line in infile:
                component = line.strip().split("|")
                if component[0] != 'mz':
                    sql = "INSERT INTO sample_component(component_id, component_name, exp_mass_charge, exp_rt, experiment_id) " \
                          "VALUES (NULL, " + "'" + component[2].strip()  + "','" + component[0].strip() + "'," + "'" + component[1].strip() + "','" + \
                          str(prot_exp_id) + "')"

                    try:
                        self.cursor.execute(sql)
                    except Exception as e:
                        self.logger.exception(e)

            self.db.commit()

        # insert hela
        sql = "INSERT INTO sample_component(component_id, component_name, exp_mass_charge, exp_rt, experiment_id) " \
              "VALUES (NULL,'Hela Digest','-1','-1','" +str(prot_exp_id) + "')"

        try:
            self.cursor.execute(sql)
        except Exception as e:
            self.logger.exception(e)

        self.db.commit()

    def insert_machines(self):
        # id, name, (serial), (description), venue, use_metab, use_prot
        # machines
        with open(self.fs.config_dir + "\\" +'machines.txt', 'r') as infile:
            for line in infile:
                machine = line.strip().split("|")

                sql = "INSERT INTO machine(machine_id, machine_name, machine_venue, use_metab, use_prot) " \
                      "VALUES (NULL, " + "'" + machine[0].strip() + "','" + machine[3].strip() + \
                      "','" + machine[4].strip() + "','" + machine[5].strip() +"')"

                try:
                    self.cursor.execute(sql)
                except Exception as e:
                    self.logger.exception(e)

            self.db.commit()

    def insert_experiments(self):

            sql = "INSERT INTO experiment(experiment_id, experiment_type) " \
                  "VALUES (NULL, " + "'proteomics')"
            try:
                self.cursor.execute(sql)
            except Exception as e:
                self.logger.exception(e)

            sql = "INSERT INTO experiment(experiment_id, experiment_type) " \
                  "VALUES (NULL, " + "'metabolomics')"
            try:
                self.cursor.execute(sql)
            except Exception as e:
                self.logger.exception(e)

            self.db.commit()

    def insert_mzmine_metrics(self):
        # name,full name,disp. order, disp. name, use_metab, use_prot
        with open(self.fs.config_dir + "\\" + 'mzmine_metrics.txt', 'r') as infile:
            for line in infile:
                in_data = line.strip().split('|')
                sql = "INSERT INTO metric VALUES (NULL," + "'" + in_data[0].strip() + "'," + "'" \
                + in_data[1].strip() + "','" + in_data[2].strip() + "','" + in_data[3].strip() + \
                      "','" + in_data[4].strip() + "','" + in_data[5].strip() + "','" + "N')"

                try:
                    self.cursor.execute(sql)
                except Exception as e:
                    self.logger.exception(e)
            self.db.commit()

    def insert_morpheus_metrics(self):
        # name,full name,disp. order, disp. name, use_metab, use_prot
        with open(self.fs.config_dir + "\\" + 'morpheus_metrics.txt', 'r') as infile:
            for line in infile:
                in_data = line.strip().split('|')
                sql = "INSERT INTO metric VALUES (NULL," + "'" + in_data[0].strip() + "'," + "'" \
                      + in_data[1].strip() + "','" + in_data[2].strip() + "','" + in_data[3].strip() + \
                      "','" + in_data[4].strip() + "','" + in_data[5].strip() + "','" + "Y')"

                try:
                    self.cursor.execute(sql)
                except Exception as e:
                    self.logger.exception(e)
            self.db.commit()

    # SELECTS
    def select_and_display_all(self, tablename):
        print(tablename.upper())
        sql = "SELECT * FROM " + tablename
        
        try:
            self.cursor.execute(sql)
            result =  self.cursor.fetchall()
            for line in result:
                print(line)
        except Exception as e:
            self.logger.exception(e)

    # GETS
    def get_run_id(self, datafile):
        sql = "SELECT run_id FROM qc_run WHERE file_name = " + "'" + datafile + "'"
        try:
            self.cursor.execute(sql)
            run_id = self.cursor.fetchone()
            return run_id[0]
        except Exception as e:
            print(e)
            return False

    def get_metric_id(self, name):
        sql = "SELECT metric_id FROM metric WHERE metric_name = " + "'" + name + "'"
        try:
            self.cursor.execute(sql)
            metric_id = self.cursor.fetchone()
            return metric_id[0]
        except Exception as e:
            print(e)
            return False

    def get_measurement(self, cid, mid, rid):
        sql = "SELECT m.metric_name, c.component_name, r.value, q.date_time FROM " \
                "metric m, measurement r, qc_run q, sample_component c WHERE " \
                "r.metric_id = m.metric_id AND " \
                "r.component_id = c.component_id AND " \
                "r.run_id = q.run_id AND " \
                "c.component_id = " + str(cid) + " AND " \
                "q.run_id = " + str(rid) + " AND " \
                "m.metric_id = " + str(mid)

        try:
            self.cursor.execute(sql)
            measurement = self.cursor.fetchone()
            return measurement
        except Exception as e:
            print(e)
            return False

if __name__ == "__main__":

    user = "root"
    password = "metabolomics"
    database = "mpmfdb"
    fs = FileSystem("", "", "", "")

    db = MPMFDBSetUp(user, password, database, fs)

    if db.connected:
        db.set_up()
        db.logger.warning("Database Loaded Successfully")







