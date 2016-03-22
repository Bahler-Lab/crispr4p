#!/usr/bin/python
# Import modules for CGI handling
import cgi, cgitb
import os, json
import crispr4p.crispr4p as crp


class viewForm(object):
    def __call__(self):
        text = self.print_html_header()
        text += self.load_bahler_template()
        text += self.generate_form()
        print text

    def print_html_header(self):
        return "Content-type: text/html\n\n"

    def load_bahler_template(self):
        src_path = os.path.dirname(__file__)
        try:
            if src_path == "":
                src_path = "."
            template_file = open(src_path + '/template/bahler_template.html')
        except IOError as err:
            print err
        return template_file.read()

    def generate_form(self):
        form_text = '<H3>Please, enter gene information.</H3>\
                <FORM METHOD = post ACTION = \"webapp.py\">\
                You can either specify by gene name,<br>\
                Gene name: <INPUT type = text name = \"name\"><br><br><br>\
                or specify the coordinates:<br>\
                Chromosome: <input type = "text" name="chromosome" >\
                Coordinates: from <INPUT type="number" name="coor_lower" min="0" max="1000000">\
                to <INPUT type = number name = "coor_upper" min="0" max="1000000"><br>\
                e.g. chromosome = I; Coordinates from 112000 to 115000<br><br>\
                <input type="submit" name="action" value="Submit"><br><br>\
                </FORM>'
        return form_text



class PrimerDesignModel(object):
    def __init__(self, name=None, cr=None, start=None, end=None):
        self.name = name
        self.cr = cr
        self.start = start
        self.end = end
    
    def check_webinput(self):
        print 'You are cheking: <br>'
        if self.name==None:
            check_txt = 'chromosome = ' + self.cr + '<br>'
            check_txt +='coordinates = from' + self.start + 'to' + self.end + '<br>'
            return check_txt
        else:
            check_txt = 'name = ' + self.name + '<br>'
            return check_txt

    def run(self):
        datapath = "data/"
        FASTA = datapath + 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
        COORDINATES = datapath + 'COORDINATES.txt'
        SYNONIMS = datapath + 'SYNONIMS.txt'

        pd = crp.PrimerDesign(FASTA, COORDINATES, SYNONIMS)
        try:
            self.tablePos_grna, self.hr_dna, self.primercheck = pd.runWeb(self.name, self.cr, self.start, self.end, nMismatch=0)
        except AssertionError as err:
            print "Error: ", err

    def result_html(self):
        result_html = self.check_webinput()

        left_table = self.gRNA_table(self.tablePos_grna)


        right_table = self.HR_DNA_report(self.hr_dna)
        right_table += self.CheckingPrimers_report(self.primercheck)

        right_table = '''<div style="float:right;width:400;">
        <table style="table-layout: fixed" ,="" border="0">
        %s%s
        </table></div>
        ''' % (self.show_Primer(), right_table)

        divider = '%s%s'
        divider = divider % (left_table, right_table)

        divider = divider + '</body></html>'

        return result_html+divider

    
    def write_html_table(self, table, name):
        html_table = '''
        <tbody><tr>
                <th align="right">''' + name + '''</th>
            </tr>'''
        for row in table:
            html_table += '<td "="" align="right" bgcolor="Azure" valign="top">' + row[0] + "</td>"
            html_table +=  '<td style="width: 80%; word-wrap: break-word" bgcolor="SeaShell">' + str(row[1]) + '</td></tr>'
        html_table += '''</tbody>'''
        return html_table

    def gRNA_table(self, table):

        table_script = '''
        <div id="leftTable" style="float:left;width:500;overflow: scroll; height:600;">
        <table id="posTable" style=" table-layout: fixed; border:0;">
            <tbody><tr>
                <th>Pos</th><th style="width:60%;">gRNA</th>
                <th>8</th><th>10</th><th>12</th><th>14</th><th>16</th><th>18</th><th>20</th>
            </tr>
        </tbody></table><br>
        <form action="">
                        <script>

                        function update_primer(number)
                        {
                            var primer = rows[number][1]
                            document.getElementById('primer_algo').innerHTML = primer[0];
                            document.getElementById('primer_forward').innerHTML = primer[1];
                            document.getElementById('primer_reverse').innerHTML = primer[2];
                            document.getElementById('primer_pos').innerHTML = primer[3];
                            document.getElementById('primer_strand').innerHTML = primer[4];
                        }

                        var rows = ''' + json.dumps(table) + ''';
                        function createTable() {
                            var table = document.getElementById("posTable");
                            var len = 0;
                            for (i = 0, len = rows.length; i < len; i++){
                                var row = table.insertRow(i+1);

                                var radioHtml = '<input type="radio" name="gRNA" value="' + i.toString() + '" onclick="update_primer(this.value);"';
                                if ( i == 0 ) {
                                    radioHtml += ' checked="checked"';
                                }
                                    radioHtml += '/>';

                                row.insertCell(0).innerHTML = radioHtml;

                                row.insertCell(1).innerHTML = rows[i][0];
                                row.insertCell(2).innerHTML = rows[i][2];
                                row.insertCell(3).innerHTML = rows[i][3];
                                row.insertCell(4).innerHTML = rows[i][4];
                                row.insertCell(5).innerHTML = rows[i][5];
                                row.insertCell(6).innerHTML = rows[i][6];
                                row.insertCell(7).innerHTML = rows[i][7];

                            }
                        }

                        window.onload = function(){createTable();update_primer(0);}
                        </script>
        </form>
        </div>
        '''


        return table_script

    def show_Primer(self):
        primer_dict = [
            ('gRNA', '<div id="primer_algo"></div>'),
            ('Forward:', '<div id="primer_forward"></div>'),
            ('Reverse:', '<div id="primer_reverse"></div>'),
            ('Pos:', '<div id="primer_pos"></div>'),
            ('Strand:', '<div id="primer_strand"></div>')
        ]
        return '''<div id="primer">%s<div id="primer" />''' % self.write_html_table(primer_dict, 'Primer')

    def gRNA_report(self, gRNA, start):
        '''
        Not in use
        :param gRNA:
        :param start:
        :return:
        '''
        gRNA_dict = [('gRNA:',   gRNA[0]),
                ('PAM: ',  "Position=" + str(int(gRNA[3])+int(start)) + \
                        " sequence=" + str(gRNA[5]) + " strand=(" + str(gRNA[4]) + ")"),
                ('gRNA-fw: ', gRNA[1]),
                ('gRNA-rv: ', gRNA[2])]
        return self.write_html_table(gRNA_dict, "gRNA")

    def HR_DNA_report(self, hr_dna):
        hr_dna_dict = [('HRfw: ', hr_dna[0]),
                ('HRrv: ', hr_dna[1]),
                ('Deleted DNA: ', hr_dna[2])]
        return self.write_html_table(hr_dna_dict, "Homologous Recombination Primers")

    def CheckingPrimers_report(self, primerDesigns):
        pm = primerDesigns[0]
        pm_dict = [
                ('Check primer left: ', pm['PRIMER_LEFT_0_SEQUENCE']), 
                ('LEFT_TM:', int(round(pm['PRIMER_LEFT_0_TM']))),
                ('Check primer right: ', pm['PRIMER_RIGHT_0_SEQUENCE']), 
                ('RIGHT_TM:', int(round(pm['PRIMER_RIGHT_0_TM']))),
                ('Deleted DNA product size: ', str(pm['PRIMER_PAIR_0_PRODUCT_SIZE']) + "(bp)"),
                ('Negative result product size: ', str(pm['negative_result']) + "(bp)")]

        return self.write_html_table(pm_dict, "Checking Primers")

class controller(object):
    def __init__(self):
        self.form = cgi.FieldStorage()
    
    def check_form_action(self):
        #todo: change this

        coordinate_form_names = ["coor_upper", "coor_lower", "chromosome"]
        get_form_val = lambda x: str(self.form.getvalue(x))

        if self.form.has_key("action"):
            if self.form.has_key("name"):
                return [get_form_val("name"), None, None, None]
            elif all(j==True for j in [self.form.has_key(i) for i in coordinate_form_names]):
                return [None, get_form_val("chromosome"), 
                        get_form_val("coor_lower"),
                        get_form_val("coor_upper")]
            else:
                print '<font color="red"> Error: neither inputs of name mode nor\
                        coordinate mode is complete</font>' # for web user
                raise ValueError("neither inputs of name mode nor coordinate\
                        mode is complete.")
        else: return None

    def is_render(self):
        if self.form.has_key("render"): return True
        else: return False


    def run_model(self):
        if not self.is_render():
            model_arguments = self.check_form_action()
            if model_arguments != None:
                self.model = PrimerDesignModel(*model_arguments)
                self.model.run()
                print self.model.result_html()


def webrun():
    init_form = viewForm()
    init_form()
    model = controller()
    model.run_model()




webrun()


