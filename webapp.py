#!/usr/bin/python
# Import modules for CGI handling
import cgi
import os, json
import crispr4p.crispr4p as crp


class viewForm(object):
    def __call__(self):
        text = self.print_html_header()
        text += self.load_bahler_template()
        return text

    def print_html_header(self):
        return "Content-type: text/html\n\n"

    def load_bahler_template(self):
        src_path = os.path.dirname(__file__)
        src_path = "." if src_path == "" else src_path

        try:
            with open(src_path + '/template/bahler_template.html') as fh:
                template_file = fh.read()
            return template_file
        except IOError as err:
            print err


class PrimerDesignModel(object):
    def __init__(self, name=None, cr=None, start=None, end=None):
        self.name = name
        self.cr = cr
        self.start = start
        self.end = end
        self.primercheck = None

    def run(self):
        datapath = "data/"
        FASTA = datapath + 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
        COORDINATES = datapath + 'COORDINATES.txt'
        SYNONIMS = datapath + 'SYNONIMS.txt'

        pd = crp.PrimerDesign(FASTA, COORDINATES, SYNONIMS, precomputed_folder='precomputed')
        try:
            self.tablePos_grna, self.hr_dna, self.primercheck, self.name, self.cr, self.start, self.end = pd.runWeb(self.name, self.cr, self.start, self.end, nMismatch=0)
        except AssertionError as err:
            print "Error: ", err

    def result_html(self):
        pm = self.primercheck[0] if self.primercheck else {}
        result_dict = {'name': self.name,
                       'chromosome': self.cr,
                       'start': self.start,
                       'end': self.end,
                       'hrfw': self.hr_dna[0],
                       'hrrv': self.hr_dna[1],
                       'deleted_dna': self.hr_dna[2],
                       'primer_left': pm.get('PRIMER_LEFT_0_SEQUENCE', '-'),
                       'left_tm': "%d &deg;C" % int(round(pm.get('PRIMER_LEFT_0_TM', '0'))),
                       'primer_right': pm.get('PRIMER_RIGHT_0_SEQUENCE', '-'),
                       'right_tm': "%d &deg;C" %  int(round(pm.get('PRIMER_RIGHT_0_TM', '0'))),
                       'deleted_dna_size': str(pm.get('PRIMER_PAIR_0_PRODUCT_SIZE', '-')) + " (bp)",
                       'negative_result_size': str(pm.get('negative_result', '-')) + " (bp)"}

        result_dict['json_table'] = json.dumps(self.tablePos_grna)

        src_path = os.path.dirname(__file__) if os.path.dirname(__file__) else '.'
        with open(src_path + '/template/container_table.html') as fh:
            template_file = fh.read()
        return template_file % (result_dict)


class controller(object):
    def __init__(self):
        self.form = cgi.FieldStorage()

    def check_form_action(self):

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
                try:
                    self.model.run()
                    ans = self.model.result_html()
                except:
                    ans = '<font color="red"><h2>ERROR: please contact to: <a href="mailto:m.rodriguezlopez@ucl.ac.uk">m.rodriguezlopez@ucl.ac.uk</a></h2></font>'
                finally:
                    return ans
        return  ''


def webrun():
    init_form = viewForm()
    temp = init_form()
    model = controller()
    ans = model.run_model()
    print temp % ans


webrun()
