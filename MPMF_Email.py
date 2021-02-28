import smtplib,ssl
from MPMF_File_System import FileSystem
from MPMF_Database_SetUp import MPMFDBSetUp
from email.message import EmailMessage
from bs4 import BeautifulSoup


class SendEmail:
    """
        Email module
        Creates and sends emails
        Used by ProcessRawFile
        https://realpython.com/documenting-python-code/#documenting-your-python-code-base-using-docstrings
    """
    def __init__(self, limit_data, db):

        self.port = 465  # For SSL
        self.smtp_server = "smtp.gmail.com"
        self.sender_email = "mpmf.qc@gmail.com"
        self.password = "Proteomics1!"
        self.context = ssl.create_default_context()
        self.email_data = limit_data
        self.db = db
        self.date_time = self.get_date_time()
        self.subject = self.create_subject()
        #self.contacts_file = "Config\contacts\contacts-" +  self.email_data['metadata']['loc'].lower() + "-" + \
                             #self.email_data['metadata']['experiment'].lower() + ".csv"
        self.contacts_file = "Config\contacts\contacts.csv"
        self.send_message()

    def get_date_time(self):
        sql = "SELECT date_time FROM qc_run WHERE file_name = " + "'" +  self.email_data['metadata']['filename'] + "'"
        self.db.cursor.execute(sql)
        dt = self.db.cursor.fetchone()[0]
        return dt

    def send_message(self):

        with smtplib.SMTP_SSL(self.smtp_server, self.port, context=self.context) as server:
            server.login(self.sender_email, self.password)
            with open(self.contacts_file) as file:
                addresses = file.readlines()
            addresses = [address.split("|") for address in addresses]


            for email in addresses:
                try:
                    msg = EmailMessage()
                    msg.set_content(self.create_body(email[1].strip()))
                    msg.add_alternative(self.create_html_body(email[1].strip()), subtype='html')
                    msg['Subject'] = self.subject
                    msg['From'] = self.sender_email
                    msg['To'] = email[0].strip()
                    server.send_message(msg)
                except Exception as e:
                    print(email[0])
                    print(e)

    def create_subject(self):
        return "WARNING:" + self.email_data['metadata']['experiment'].upper() + " QC for " + \
               self.email_data['metadata']['machine']

    def create_body(self, name):
        intro = "Hi " + name + ",\n\n" + \
            "This is the automated " + self.email_data['metadata']['experiment'] + " QC " + \
            "for " + self.email_data['metadata']['machine'] + ".\n\n" + \
            "Insufficient counts were encounterd on \n" + str(self.date_time) + "\n\n" + \
            "For the following metrics:\n\n"

        message = ""
        for metric in self.email_data:
            if metric != 'metadata':
                # set display metric names
                temp = metric
                if metric == 'tf':
                    metric = 'TAILING FACTOR'
                elif metric == 'af':
                    metric = 'ASYMMETRY FACTOR'
                elif metric == 'fwhm':
                    metric = 'FULL WIDTH HALF MAXIMUM'
                elif metric == 'rt':
                    metric = 'RETENTION TIME'
                message += metric.upper() + "\n"
                metric = temp
                for comp in self.email_data[metric]:
                    message += comp + "\t" + str(self.email_data[metric][comp]) + "\n"
            message += "\n"

        return intro + message

    def create_html_body(self, name):
        with open("Config\email_template.html") as fp:
            self.soup = BeautifulSoup(fp,features="html.parser")

        # fill in intro
        tag_name = self.soup.find(id="name")
        tag_name.string = name
        tag_exp = self.soup.find(id="exp")
        tag_exp.string = self.email_data['metadata']['experiment'].upper()
        tag_mach = self.soup.find(id="machine")
        tag_mach.string = self.email_data['metadata']['machine'].upper()
        tag_dt = self.soup.find(id="datetime")
        tag_dt.string = str(self.date_time.strftime("%A, %d %B %Y %I:%M%p"))

        # add tables
        tag_tables = self.soup.find(id="tables")
        for metric in self.email_data:
            if metric != 'metadata':
                # set display metric names
                temp = metric
                if metric == 'tf':
                    metric = 'TAILING FACTOR'
                elif metric == 'af':
                    metric = 'ASYMMETRY FACTOR'
                elif metric == 'fwhm':
                    metric = 'FULL WIDTH HALF MAXIMUM'
                elif metric == 'rt':
                    metric = 'RETENTION TIME'

                # table header
                new_tag = self.soup.new_tag("h3")
                bold_tag = self.soup.new_tag("strong")
                bold_tag.string = metric.upper()
                new_tag.append(bold_tag)
                tag_tables.append(new_tag)

                metric = temp
				
				# table content
                table_tag = self.soup.new_tag("table")
                table_tag['style'] = "background-color:blanchedalmond"
                for comp in self.email_data[metric]:
                    row_tag = self.soup.new_tag("tr")
                    comp_tag = self.soup.new_tag("td")
                    i_tag = self.soup.new_tag("strong")
                    i_tag.string = comp
                    comp_tag.append(i_tag)
                    row_tag.append(comp_tag)
                    if self.email_data[metric][comp][0] != "NO VALUE":
                        for entry in range(len(self.email_data[metric][comp])):
                            col_tag = self.soup.new_tag("td")
                            col_tag.string = self.email_data[metric][comp][entry]
                            row_tag.append(col_tag)
                        table_tag.append(row_tag)
                    else:
                        col_tag = self.soup.new_tag("td")
                        col_tag['colspan'] = 2
                        col_tag.string = "NO VALUE"
                        row_tag.append(col_tag)
                        table_tag.append(row_tag)
                tag_tables.append(table_tag)
                #print(self.email_data)

        #print(self.soup.prettify())
        return str(self.soup)


# include functions for testing
def check_email_thresholds_prot(filename):
        # checks metric values against the thresholds in config files
        # and sends email if any outsdide limits

        # get threshold limits
        with open(fs.thresh_email) as f:
            limits = f.readlines()

        # remove header
        limits.pop(0)

        # create dict for storing
        thresholds = {}
        for limit in limits:
            new_limit = limit.split("|")
            if new_limit[1] != '':
                thresholds[new_limit[0]] = [new_limit[1], new_limit[2], new_limit[3].strip()]

        # get run_id
        sql = "SELECT * FROM qc_run WHERE file_name =" + "'" + filename + "'"
        db.cursor.execute(sql)
        run_id = db.cursor.fetchone()[0]

        breaches = {}
        for metric in thresholds:

            # limits
            tot = int(thresholds[metric][0])
            lower = thresholds[metric][1]
            upper = thresholds[metric][2]

            # get metric_id
            sql = "SELECT metric_id FROM metric WHERE metric_name = " + "'" + metric + "'"
            db.cursor.execute(sql)
            metric_id = db.cursor.fetchone()[0]

            # get values for metric and run_id
            sql = "SELECT c.component_name, v.value FROM " + \
                  "measurement v, sample_component c, metric m " + \
                  "WHERE m.metric_id = v.metric_id AND " + \
                  "c.component_id = v.component_id AND " + \
                  "v.run_id = " + "'" + str(run_id) + "'" + \
                  " AND v.metric_id = " + "'" + str(metric_id) + "'"

            db.cursor.execute(sql)
            results = db.cursor.fetchall()

            # get components that exceed limits for each metric
            comps = {}
            if metric == "mass_error_ppm":
                # check limits
                for result in results:
                    if result[1] > float(upper) or result[1] < float(lower):
                        comps[result[0]] = [str(round(result[1], 3)) + " ppm"]
                    # catch missed values
                    if result[1] == -1000000.0:
                        comps[result[0]] = ["NO VALUE"]
            elif metric == "area_normalised":
                # check limits
                for result in results:
                    if result[1] < float(lower):
                        comps[result[0]] = [str(round(result[1], 3))]
                    if result[1] == -100.0:
                        comps[result[0]] = ["NO VALUE"]
            elif metric == "fwhm":
                # check limits
                for result in results:
                    if result[1] > float(upper):
                        comps[result[0]] = [str(round(result[1], 2)) + " sec"]
                    if result[1] == 0:
                        comps[result[0]] = ["NO VALUE"]
            elif metric == "tf":
                # check limits
                for result in results:
                    if result[1] > float(upper):
                        comps[result[0]] = [str(round(result[1], 3))]
                    if result[1] == 0:
                        comps[result[0]] = ["NO VALUE"]
            elif metric == "af":
                # check limits
                for result in results:
                    if result[1] > float(upper):
                        comps[result[0]] = [str(round(result[1], 3))]
                    if result[1] == 0:
                        comps[result[0]] = ["NO VALUE"]
            elif metric == "MS/MS Spectra":
                # determine percentiles
                sql = "SELECT value FROM measurement WHERE metric_id = " + "'" + str(metric_id) + "'" + \
                      " ORDER by value"
                db.cursor.execute(sql)
                all_results = db.cursor.fetchall()
                all_values = [float(item[0]) for item in all_results]

                # get index in ordered list of values
                try:
                    pos = all_values.index(float(results[0][1]))
                    # check upper percentile
                    if (1 - pos / len(all_values)) < float(upper) / 100:
                        comps[metric] = [str(int(results[0][1])),"Top " + str(round((1-pos / len(all_values))*100, 2)) + "%"]
                    # check lower percentile
                    elif pos / len(all_values) < float(lower) / 100:
                        comps[metric] = [str(int(results[0][1])),"Bottom " + str(round((pos / len(all_values))*100, 2)) + "%"]
                except ValueError:
                    pass
            elif metric == "Target PSMs":
                # determine percentiles
                sql = "SELECT value FROM measurement WHERE metric_id = " + "'" + str(metric_id) + "'" + \
                      " ORDER by value"
                db.cursor.execute(sql)
                all_results = db.cursor.fetchall()
                all_values = [float(item[0]) for item in all_results]

                # get index in ordered list of values
                try:
                    pos = all_values.index(float(results[0][1]))
                    # check lower percentile
                    if pos / len(all_values) < float(lower) / 100:
                        comps[metric] = [str(int(results[0][1])),"Bottom " + str(round((pos / len(all_values))*100, 2)) + "%"]
                except ValueError:
                    pass
            elif metric == "Unique Target Peptides":
                # determine percentiles
                sql = "SELECT value FROM measurement WHERE metric_id = " + "'" + str(metric_id) + "'" + \
                      " ORDER by value"
                db.cursor.execute(sql)
                all_results = db.cursor.fetchall()
                all_values = [float(item[0]) for item in all_results]

                # get index in ordered list of values
                try:
                    pos = all_values.index(float(results[0][1]))
                    # check lower percentile
                    if pos / len(all_values) < float(lower) / 100:
                        comps[metric] = [str(int(results[0][1])),"Bottom " + str(round((pos / len(all_values))*100, 2)) + "%"]
                except ValueError:
                    pass
            elif metric == "Target Protein Groups":
                # determine percentiles
                sql = "SELECT value FROM measurement WHERE metric_id = " + "'" + str(metric_id) + "'" + \
                      " ORDER by value"
                db.cursor.execute(sql)
                all_results = db.cursor.fetchall()
                all_values = [float(item[0]) for item in all_results]

                # get index in ordered list of values
                try:
                    pos = all_values.index(float(results[0][1]))
                    # check lower percentile
                    if pos / len(all_values) < float(lower) / 100:
                        comps[metric] = [str(int(results[0][1])),"Bottom " + str(round((pos / len(all_values))*100, 2)) + "%"]
                except ValueError:
                    pass
            elif metric == "Precursor Mass Error":
                # check limits
                for result in results:
                    if result[1] > float(upper) or result[1] < float(lower):
                        comps[metric] = [str(round(result[1], 3)) + " ppm"]
            elif metric == 'rt':
                for result in results:
                    sql = "SELECT component_id FROM sample_component WHERE component_name = " + "'" + str(
                        result[0]) + "'"
                    db.cursor.execute(sql)
                    comp_id = db.cursor.fetchone()[0]

                    # get all values per component
                    sql = "SELECT value FROM measurement WHERE metric_id = " + "'" + str(metric_id) + "'" + \
                          " AND component_id = " + "'" + str(comp_id) + "'" + " AND value <> 0 " \
                                                                              " ORDER BY value"
                    db.cursor.execute(sql)
                    all_results = db.cursor.fetchall()
                    all_values = [float(item[0]) for item in all_results]

                    # get index in ordered list of values
                    try:
                        pos = all_values.index(float(result[1]))
                        # check upper percentile
                        if (1 - pos / len(all_values)) < float(upper) / 100:
                            comps[result[0]] = [str(round(result[1], 2)) +" min","Top " + str(
                                round((1 - pos / len(all_values)) * 100, 2)) + "%"]
                        # check lower percentile
                        elif pos / len(all_values) < float(lower) / 100:
                            comps[result[0]] = [str(round(result[1], 2)) + " min" ,"Bottom " + \
                                               str(round((pos / len(all_values)) * 100, 2)) + "%"]
                    except ValueError:
                        pass

                    # catch missed values
                    if result[1] == 0:
                        comps[result[0]] = ["NO VALUE"]


            # add to breaches if tot or more
            if len(comps) >= tot:
                breaches[metric] = comps

        return breaches


def check_email_thresholds_metab(filename):

        # get threshold limits
        with open(fs.thresh_email) as f:
            limits = f.readlines()

        # remove header
        limits.pop(0)

        # create dict for storing
        thresholds = {}
        for limit in limits:
            new_limit = limit.split("|")
            if new_limit[1] != '':
                thresholds[new_limit[0]] = [new_limit[1], new_limit[2], new_limit[3].strip()]

        # get run_id
        sql = "SELECT * FROM qc_run WHERE file_name =" + "'" + filename + "'"
        db.cursor.execute(sql)
        run_id = db.cursor.fetchone()[0]

        breaches = {}
        for metric in thresholds:
            # limits
            tot = int(thresholds[metric][0])
            lower = thresholds[metric][1]
            upper = thresholds[metric][2]

            # get metric_id
            sql = "SELECT metric_id FROM metric WHERE metric_name = " + "'" + metric + "'"
            db.cursor.execute(sql)
            metric_id = db.cursor.fetchone()[0]

            # get values for metric and run_id (not limited by polarity)
            sql = "SELECT c.component_name, v.value FROM " + \
                  "measurement v, sample_component c, metric m " + \
                  "WHERE m.metric_id = v.metric_id AND " + \
                  "c.component_id = v.component_id AND " + \
                  "v.run_id = " + "'" + str(run_id) + "'" + \
                  " AND v.metric_id = " + "'" + str(metric_id) + "'" + \
                  " AND c.component_name <> 'Glucose'"

            db.cursor.execute(sql)
            results = db.cursor.fetchall()

            # get components that exceed limits for each metric
            comps = {}
            if metric == "mass_error_ppm":
                modes = ['N', 'P']
                for mode in modes:
                    # get values by polarity
                    comps = {}
                    sql = "SELECT c.component_name, v.value FROM " + \
                          "measurement v, sample_component c, metric m " + \
                          "WHERE m.metric_id = v.metric_id AND " + \
                          "c.component_id = v.component_id AND " + \
                          "v.run_id = " + "'" + str(run_id) + "'" + \
                          " AND v.metric_id = " + "'" + str(metric_id) + "'" + \
                          " AND c.component_mode =" + "'" + mode + "'" + \
                          " AND c.component_name <> 'Glucose'"

                    db.cursor.execute(sql)
                    results = db.cursor.fetchall()

                    # check limits
                    for result in results:
                        if result[1] > float(upper) or result[1] < float(lower):
                            comps[result[0]] = [str(round(result[1], 3)) + " ppm"]

                        # catch missed values
                        if result[1] == -1000000.0:
                            comps[result[0]] = ["NO VALUE"]

                    # add to breaches if tot or more
                    if len(comps) >= tot:
                        if mode == 'N':
                            breaches[metric + "_Neg"] = comps
                        else:
                            breaches[metric + "_Pos"] = comps
            elif metric == 'rt':
                for result in results:
                    sql = "SELECT component_id FROM sample_component WHERE component_name = " + "'" + str(
                        result[0]) + "'"
                    db.cursor.execute(sql)
                    comp_id = db.cursor.fetchone()[0]

                    # get all values per component
                    sql = "SELECT value FROM measurement WHERE metric_id = " + "'" + str(metric_id) + "'" + \
                          " AND component_id = " + "'" + str(comp_id) + "'" + " AND value <> 0 " \
                                                                              " ORDER BY value"
                    db.cursor.execute(sql)
                    all_results = db.cursor.fetchall()
                    all_values = [float(item[0]) for item in all_results]

                    # get index in ordered list of values
                    try:
                        pos = all_values.index(float(result[1]))
                        # check upper percentile
                        if (1 - pos / len(all_values)) < float(upper) / 100:
                            comps[result[0]] = [str(round(result[1], 2)) + " min", "Top " + str(round((1-pos / len(all_values))*100, 2)) + "%"]
                        # check lower percentile
                        elif pos / len(all_values) < float(lower) / 100:
                            comps[result[0]] = [str(round(result[1], 2)) + " min", "Bottom " + \
                                               str(round((pos / len(all_values))*100, 2))  + "%"]
                    except ValueError:
                        pass

                    # catch missed values
                    if result[1] == 0:
                        comps[result[0]] =["NO VALUE"]

                if len(comps) >= tot:
                    breaches[metric] = comps
            elif metric == 'area_normalised':
                for result in results:
                    sql = "SELECT component_id FROM sample_component WHERE component_name = " + "'" + str(
                        result[0]) + "'"
                    db.cursor.execute(sql)
                    comp_id = db.cursor.fetchone()[0]

                    # get all values per component
                    sql = "SELECT value FROM measurement WHERE metric_id = " + "'" + str(metric_id) + "'" + \
                          " AND component_id = " + "'" + str(comp_id) + "'" + " AND value <> -100 " \
                                                                              " ORDER BY value"
                    db.cursor.execute(sql)
                    all_results = db.cursor.fetchall()
                    all_values = [float(item[0]) for item in all_results]

                    # get index in ordered list of values
                    try:
                        pos = all_values.index(float(result[1]))
                        # check upper percentile
                        if (1 - pos / len(all_values)) < float(upper) / 100:
                            comps[result[0]] = [str(round(result[1], 2)) , "Top " + str(round((1-pos / len(all_values))*100, 2)) + "%"]
                        # check lower percentile
                        elif pos / len(all_values) < float(lower) / 100:
                            comps[result[0]] = [str(round(result[1], 2)) , "Bottom " + str(round((pos / len(all_values))*100, 2)) + "%"]
                    except ValueError:
                        pass

                    # catch missed values
                    if result[1] == -100:
                        comps[result[0]] = ["NO VALUE"]

                if len(comps) >= tot:
                    breaches[metric] = comps

        return breaches


if __name__ == "__main__":
    # TESTING ONLY
    exp = "metabolomics"
    e_id = '2'
    machine = "qeclassic"
    venue = "clayton"

    db_info = {"user": "root", "password": "metabolomics", "database": "mpmfdb"}
    fs = FileSystem("", "*", "", exp)
    db = MPMFDBSetUp(db_info["user"], db_info["password"], db_info["database"], fs)

    sql = "SELECT file_name FROM qc_run WHERE machine_id = (SELECT machine_id FROM machine WHERE machine_name = " + \
            "'" + machine + "'" + ") AND experiment_id = " +  "'" + e_id + "'"

    db.cursor.execute(sql)
    all_files = db.cursor.fetchall()
    lim = 0
    for filename in all_files:
        metadata = {'filename':filename[0], 'experiment':exp, 'machine':machine, 'loc': venue}
        if exp == 'metabolomics':
            email_data = check_email_thresholds_metab(filename[0])
        else:
            email_data = check_email_thresholds_prot(filename[0])
            #print(email_data)
        if len(email_data) > 0:
            email_data['metadata'] = metadata
            SendEmail(email_data, db)
        else:
            print("NO BREACHES")
        if lim > 5:
            break
        lim += 1

# NOTE: Add smpt_server, address and password to a config file
# NOTE: less secure apps needs to be set