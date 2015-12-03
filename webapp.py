#!/usr/local/bin/python2.7
# Import modules for CGI handling
import cgi, cgitb
import os
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
                <INPUT type="checkbox" name=allgene> show all gRNAs &nbsp&nbsp&nbsp\
                <input type="submit" name="action" value="Submit"><br><br>\
                </FORM>'
        return form_text



class PrimerDesignModel(object):
    def __init__(self, name=None, cr=None, start=None, end=None, allGRNA=None):
        self.name = name
        self.cr = cr
        self.start = start
        self.end = end
        self.allGRNA = allGRNA
    
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
            self.ansTuple = pd.runWeb(self.name, self.cr, self.start, self.end, allGRNA=self.allGRNA)
        except AssertionError as err:
            print "Error: ", err



    def result_html(self):
        result_html = self.check_webinput()
        if self.allGRNA:
            result_html += self.all_gRNA_report(self.ansTuple)
        else:
            result_html += self.gRNA_report(self.ansTuple[0], self.ansTuple[-1][1])
            result_html += self.HR_DNA_report(self.ansTuple[1])
            result_html += self.CheckingPrimers_report(self.ansTuple[2])
        return result_html

    
    def write_html_table(self, table, name):
        html_table = '<table style="width: 100%; table-layout: fixed", boarder="0"><tr><th align="left">' + name + '</th></tr>'
        for row in table:
            html_table += '<tr><td valign="top" align="right" bgcolor="Azure"">' + row[0] + "</td>"
            html_table +=  '<td  style="width: 80%; word-wrap: break-word" bgcolor="SeaShell">' + str(row[1]) + '</td></tr>'
        html_table += '</table><br>'
        return html_table


    def check_box_js(self):
        html = """
            <script>
            function handleClick(text, checkbox){
                if(checkbox.checked==true)
                    text.style.color='red';
                else
                    text.style.color='black'
            }
            </script>
            """
        return html

    def all_gRNA_report(self, data):
        report_html = self.check_box_js()
        count = 1
        for d in data:
            value = 'GRNA:%s; PAM:%s; Strand:%s; Start:%s;'%(d[0], d[5], d[4], d[3])
            report_html +=  '<INPUT type="checkbox" onClick="handleClick(gene%s, this)">&nbsp'%str(count)
            report_html += '<input type="text" id=gene%s value="%s" style="border:none; border-color:transparent" \
                            readonly size=100>'%(str(count), value)
            report_html += '<br>'
            count += 1
        return report_html
    
    def gRNA_report(self, gRNA, start):
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
        if self.form.has_key("allgene"): self.allgene = True
        else: self.allgene = False

        coordinate_form_names = ["coor_upper", "coor_lower", "chromosome"]
        get_form_val = lambda x: str(self.form.getvalue(x))

        if self.form.has_key("action"):
            if self.form.has_key("name"):
                return [get_form_val("name"), None, None, None, self.allgene]
            elif all(j==True for j in [self.form.has_key(i) for i in coordinate_form_names]):
                return [None, get_form_val("chromosome"), 
                        get_form_val("coor_lower"),
                        get_form_val("coor_upper"), self.allgene]
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

