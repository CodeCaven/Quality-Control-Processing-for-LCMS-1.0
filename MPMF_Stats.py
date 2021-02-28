from MPMF_File_System import FileSystem
from MPMF_Database_SetUp import MPMFDBSetUp
import pandas as pd
import numpy as np

class Stat:
    """
        Creates stats for history plots
        and the normalised metrics
        Inserts these into database and updates processing status (completed='Y')
        Used by ProcessRawFile after a batch of raw files are processed
        https://realpython.com/documenting-python-code/#documenting-your-python-code-base-using-docstrings
    """
    def __init__(self, exp_type, db, machine):
        self.e_type = exp_type
        self.db = db
        self.machine = machine

    def run(self):
        if self.check_for_updates():
            self.compute_stats()
            self.insert_update_stats()
            self.compute_derived_metrics()
            print("Updated stats and normalised metrics for mzmine " + self.machine + " " + self.e_type)
            if self.e_type == 'PROTEOMICS':
                self.compute_morpheus_stats()
                self.insert_update_morpheus_stats()
                print("Updated morpheus stats " + self.machine + " " + self.e_type)
            self.update_completed()

        else:
            print("No new data for " + self.machine + " " + self.e_type)

    def check_for_updates(self):

        sql = "SELECT COUNT(*) FROM qc_run WHERE machine_id = " + \
                " (SELECT machine_id FROM machine WHERE machine_name = '" + \
                self.machine + "') AND experiment_id = (SELECT experiment_id FROM experiment " + \
                " WHERE experiment_type = '" + self.e_type.lower() + "')" + " AND completed = 'N'"

        count = 0
        try:
            self.db.db.commit() # not picking up commits in processing from other processes?
            self.db.cursor.execute(sql)
            count = int(self.db.cursor.fetchone()[0])
            #print(count)
        except Exception as e:
            print(e)

        if count > 0:
            return True
        else:
            return False

    def compute_stats(self):
        # creates a pandas table of stats for the stats table

        # exclude area and height normalised...PLUS maxQuant metrics??
        if self.e_type == "METABOLOMICS":
            sql_m = "SELECT metric_name FROM metric WHERE metric_name <> '" + "area_normalised' " + \
                    "AND metric_name <> '" + "height_normalised'" + " AND use_metab = 'Y'"
        elif self.e_type == "PROTEOMICS":
            sql_m = "SELECT metric_name FROM metric WHERE metric_name <> '" + "area_normalised' " + \
                    "AND metric_name <> '" + "height_normalised'" + " AND use_prot = 'Y'" + " AND morpheus = 'N'" # limiting morpheus?

        sql_c = "SELECT component_name FROM sample_component WHERE component_name <> '" + "Glucose'" +  " AND component_name <> 'Hela Digest'" + \
                " AND experiment_id = (SELECT experiment_id FROM experiment WHERE experiment_type = '" + self.e_type.lower() + "')"
        # limiting hela?

        try:
            self.db.cursor.execute(sql_m)
            metrics = self.db.cursor.fetchall()
        except Exception as e:
            print(e)

        try:
            self.db.cursor.execute(sql_c)
            components = self.db.cursor.fetchall()
        except Exception as e:
            print(e)

        all_stats = []
        for metric in metrics:
            for component in components:
                stats = []
                sql_value = "SELECT r.value, q.date_time FROM measurement r, qc_run q, " \
                            "metric m, sample_component c, machine t " \
                            "WHERE r.metric_id = m.metric_id " \
                            "AND r.component_id = c.component_id " \
                            "AND r.run_id = q.run_id " \
                            "AND t.machine_id = q.machine_id " \
                            "AND c.component_name = " + "'" + str(component[0]) + "' " + \
                            "AND m.metric_name = " + "'" + str(metric[0]) + "' " + \
                            "AND t.machine_name = " + "'" + str(self.machine) + "'"

                try:
                    self.db.cursor.execute(sql_value)
                    values = self.db.cursor.fetchall()
                    new_line = [component[0], metric[0], self.machine]
                    for value in values:
                        # don't include missed values in stats
                        if not (float(value[0]) == 0 or float(value[0]) < -1000):
                            stats.append(float(value[0]))
                    ds = pd.Series(stats)
                    summary = ds.describe()
                    for stat in summary:
                        new_line.append(stat)
                    all_stats.append(new_line)
                except Exception as e:
                    print(e)

        # a pandas table
        self.df = pd.DataFrame(all_stats,
                          columns=['Component', 'Metric', 'Machine', 'count', 'mean', 'std', 'min', '25%',
                                   '50%', '75%', 'max'])
        #df = df.sort_values(by=['Metric', 'Venue'])
        # print(df)
        #print(len(metrics) * len(components) * len(venues))
        #df.to_csv('stats.csv')

    def compute_morpheus_stats(self):

        sql = "SELECT metric_name FROM metric WHERE use_prot ='Y' AND morpheus = 'Y'"

        try:
            self.db.cursor.execute(sql)
            metrics = self.db.cursor.fetchall()
        except Exception as e:
            print(e)

        all_stats = []
        for metric in metrics:
            stats = []
            sql_value = "SELECT r.value, q.date_time FROM measurement r, qc_run q, " \
                        "metric m, sample_component c, machine t " \
                        "WHERE r.metric_id = m.metric_id " \
                        "AND r.component_id = c.component_id " \
                        "AND r.run_id = q.run_id " \
                        "AND t.machine_id = q.machine_id " \
                        "AND c.component_name = " + "'Hela Digest'" + \
                        "AND m.metric_name = " + "'" + str(metric[0]) + "' " + \
                        "AND t.machine_name = " + "'" + str(self.machine) + "'"

            try:
                self.db.cursor.execute(sql_value)
                values = self.db.cursor.fetchall()
                new_line = ["Hela Digest", metric[0], self.machine]
                for value in values:
                    # don't include missed values in stats
                    if not (float(value[0]) == 0 or float(value[0]) < -1000):
                        stats.append(float(value[0]))
                ds = pd.Series(stats)
                summary = ds.describe()
                for stat in summary:
                    new_line.append(stat)
                all_stats.append(new_line)
            except Exception as e:
                print(e)

        # a pandas table
        self.hela_df = pd.DataFrame(all_stats,
                               columns=['Component', 'Metric', 'Machine', 'count', 'mean', 'std', 'min', '25%',
                                        '50%', '75%', 'max'])

    def insert_update_morpheus_stats(self):


        sql_c = "SELECT component_id FROM sample_component WHERE component_name = 'Hela Digest'"

        try:
            self.db.cursor.execute(sql_c)
            cid = self.db.cursor.fetchone()
        except Exception as e:
            print(e)

        sql_mac = "SELECT machine_id FROM machine WHERE machine_name = '" + self.machine + "'"

        try:
            self.db.cursor.execute(sql_mac)
            macid = self.db.cursor.fetchone()
        except Exception as e:
            print(e)

        for i in range(len(self.hela_df)):
            sql_m = "SELECT metric_id FROM metric WHERE metric_name = " + "'" + self.hela_df.iloc[i]['Metric'] + "'"

            try:
                self.db.cursor.execute(sql_m)
                mid = self.db.cursor.fetchone()
            except Exception as e:
                print(e)

            if not self.is_inserted_stat(mid[0], cid[0], macid[0]):
                sql = "INSERT into stat VALUES(" + "'" + str(mid[0]) + "','" + str(cid[0]) + "','" + \
                             str(macid[0]) + "','"  + str(self.hela_df.iloc[i]['count']) + "','" + str(
                    self.hela_df.iloc[i]['mean']) + \
                             "','" + str(self.hela_df.iloc[i]['std']) + "','" + str(self.hela_df.iloc[i]['min']) + "','" + \
                             str(self.hela_df.iloc[i]['25%']) + "','" + str(self.hela_df.iloc[i]['50%']) + "','" + str(
                    self.hela_df.iloc[i]['75%']) + \
                             "','" + str(self.hela_df.iloc[i]['max']) + "')"
            else:
                # test in next update
                sql = "UPDATE stat SET count = '" + str(self.hela_df.iloc[i]['count']) + "', mean = '" + str(self.hela_df.iloc[i]['mean']) + \
                        "', std = '" + str(self.hela_df.iloc[i]['std']) + "', min = '" + str(self.hela_df.iloc[i]['min']) + \
                        "', 25_percent = '" + str(self.hela_df.iloc[i]['25%'])  + "', 50_percent = '" + str(self.hela_df.iloc[i]['50%'])  + \
                        "', 75_percent = '" + str(self.hela_df.iloc[i]['75%']) + "', max = '" + str(self.hela_df.iloc[i]['max'])  + \
                        "' WHERE metric_id = '" + str(mid[0]) + "' AND component_id = '" + str(cid[0]) + \
                        "' AND machine_id = '" + str(macid[0]) + "'"
            try:
                self.db.cursor.execute(sql)
            except Exception as e:
                print(e)

            self.db.db.commit()

    def insert_update_stats(self):
        # run after compute_stats()
        # columns = ['Component', 'Metric', 'Machine', 'count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']

        sql_mac = "SELECT machine_id FROM machine WHERE machine_name = '" + self.machine + "'"

        try:
            self.db.cursor.execute(sql_mac)
            macid = self.db.cursor.fetchone()
        except Exception as e:
            print(e)

        for i in range(len(self.df)):
            sql_c = "SELECT component_id FROM sample_component WHERE component_name = " + \
                    "'" + self.df.iloc[i]['Component'] + "'"
            sql_m = "SELECT metric_id FROM metric WHERE metric_name = " + "'" + self.df.iloc[i]['Metric'] + "'"


            try:
                self.db.cursor.execute(sql_c)
                cid = self.db.cursor.fetchone()
            except Exception as e:
                print(e)

            try:
                self.db.cursor.execute(sql_m)
                mid = self.db.cursor.fetchone()
            except Exception as e:
                print(e)

            if not self.is_inserted_stat(mid[0], cid[0], macid[0]):
                sql = "INSERT into stat VALUES(" + "'" + str(mid[0]) + "','" + str(cid[0]) + "','" + \
                             str(macid[0]) + "','"  + str(self.df.iloc[i]['count']) + "','" + str(
                    self.df.iloc[i]['mean']) + \
                             "','" + str(self.df.iloc[i]['std']) + "','" + str(self.df.iloc[i]['min']) + "','" + \
                             str(self.df.iloc[i]['25%']) + "','" + str(self.df.iloc[i]['50%']) + "','" + str(
                    self.df.iloc[i]['75%']) + \
                             "','" + str(self.df.iloc[i]['max']) + "')"
            else:
                sql = "UPDATE stat SET count = '" + str(self.df.iloc[i]['count']) + "', mean = '" + str(self.df.iloc[i]['mean']) + \
                        "', std = '" + str(self.df.iloc[i]['std']) + "', min = '" + str(self.df.iloc[i]['min']) + \
                        "', 25_percent = '" + str(self.df.iloc[i]['25%'])  + "', 50_percent = '" + str(self.df.iloc[i]['50%'])  + \
                        "', 75_percent = '" + str(self.df.iloc[i]['75%']) + "', max = '" + str(self.df.iloc[i]['max'])  + \
                        "' WHERE metric_id = '" + str(mid[0]) + "' AND component_id = '" + str(cid[0]) + \
                        "' AND machine_id = '" + str(macid[0]) + "'"
            try:
                self.db.cursor.execute(sql)
            except Exception as e:
                print(e)

            self.db.db.commit()

    def is_inserted_stat(self, mid, cid, macid):

        sql = "SELECT * FROM stat WHERE metric_id = '" + str(mid) + "' AND " \
                "component_id = '"  + str(cid) + "' AND " \
                 "machine_id = '" + str(macid) + "'"

        try:
            self.db.cursor.execute(sql)
            result = self.db.cursor.fetchall()
            if len(result) == 0:
                return False
            else:
                return True
        except Exception as e:
           print(e)

    def compute_derived_metrics(self):
        # call after insert_update_stats

        # includes missed values to handle in Google charts

        metrics = ["area", "height"]
        for metric in metrics:
            norm_max = 0
            norm_min = 0 # test this query, limit to experiment type...join experiment and run ??
            sql = "SELECT m.metric_id, c.component_id, n.machine_id, r.run_id, r.value FROM metric m, " + \
                  "measurement r, sample_component c, qc_run q , machine n, experiment e WHERE " + \
                  "r.metric_id = m.metric_id AND " + \
                  "e.experiment_id = c.experiment_id AND " + \
                  "n.machine_id = q.machine_id AND " + \
                  "r.run_id = q.run_id AND " + \
                  "r.component_id = c.component_id AND " + \
                  "c.component_name <> '" + "Glucose' AND " + \
                  "c.component_name <> '" + "Hela Digest' AND " + \
                  "m.metric_name = '" + metric + "' AND " + \
                    "n.machine_name = '" + str(self.machine) + "' AND " + \
                    "e.experiment_type = '" + self.e_type + "'"


            try:
                self.db.cursor.execute(sql)
                values = self.db.cursor.fetchall()
                for value in values:
                    # REFACTOR: cache medians to avoid repeated db calls
                    get_median = "SELECT 50_percent FROM stat WHERE metric_id = '" + str(value[0]) + "' AND " + \
                                 "component_id = '" + str(value[1]) + "' AND " + \
                                 "machine_id = '" + str(value[2]) + "'"

                    try:
                        self.db.cursor.execute(get_median)
                        median = self.db.cursor.fetchone()
                    except Exception as e:
                        print(e)

                    # compute normalised median
                    # handle missed values by setting to -100
                    if median[0] != 0 and value[4] != 0:
                        norm = np.log2(float(value[4] / median[0]))
                    else:
                        norm = -100

                    if norm == -np.inf:
                        norm = -100

                    if norm != -100:
                        if norm > norm_max:
                            norm_max = norm
                        elif norm < norm_min:
                            norm_min = norm

                    # get metric id
                    get_mid = "SELECT metric_id FROM metric WHERE metric_name = " + "'" + metric + "_normalised" + "'"
                    try:
                        self.db.cursor.execute(get_mid)
                        mid = self.db.cursor.fetchone()
                    except Exception as e:
                        print(e)

                    # deleting means all values are updated using current stats
                    if self.is_inserted_measurement(mid[0], value[1], value[3]):
                        self.delete_measurement(mid[0], value[1], value[3])

                    insert_med = "INSERT INTO measurement VALUES('" + str(mid[0]) + "',' " + \
                                 str(value[1]) + "','" + str(value[3]) + "','" + str(norm) + "')"
                    try:
                        self.db.cursor.execute(insert_med)
                    except Exception as e:
                        print(e)
                        # print(norm)

            except Exception as e:
                print(e)

            #print("Max Norm " + metric + " = " + str(norm_max))
            #print("Min Norm " + metric + " = " + str(norm_min))

        self.db.db.commit()

    def is_inserted_measurement(self, mid, cid, rid):

        sql = "SELECT * FROM measurement WHERE metric_id = '" + str(mid) + "' AND " \
                "component_id = '"  + str(cid) + "' AND " \
                 "run_id = '" + str(rid) + "'"

        try:
            self.db.cursor.execute(sql)
            result = self.db.cursor.fetchall()
            if len(result) == 0:
                return False
            else:
                return True
        except Exception as e:
            print(e)

    def delete_measurement(self, mid, cid, rid):

        sql = "DELETE FROM measurement WHERE metric_id = '" + str(mid) + "' AND " \
                "component_id = '"  + str(cid) + "' AND " \
                 "run_id = '" + str(rid) + "'"

        try:
            self.db.cursor.execute(sql)
        except Exception as e:
            print(e)

        self.db.db.commit()

    def update_completed(self):
        sql = "UPDATE qc_run SET completed = 'Y' WHERE machine_id = " + \
            "(SELECT machine_id FROM machine WHERE machine_name = '" + \
            self.machine + "') AND experiment_id = (SELECT experiment_id FROM experiment " + \
            " WHERE experiment_type = '" + self.e_type.lower() + "')"

        try:
            self.db.cursor.execute(sql)
            self.db.db.commit()
            return True
        except Exception as e:
            print(e)


if __name__ == "__main__":

    # TESTING

    # set database details
    db_info = {"user" :"root", "password" :"metabolomics", "database" :"mpmfdb"}

    # get file system and database objects
    fs = FileSystem("Database", "", "qehfmips", "proteomics")
    db = MPMFDBSetUp(db_info["user"], db_info["password"], db_info["database"], fs)

    new_stat = Stat('PROTEOMICS', db, 'qehfmips')
    new_stat.run()
    #new_stat.run()


