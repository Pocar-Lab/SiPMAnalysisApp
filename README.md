# SiPMAnalysisApp

This is the LXe Applet and accompanying data. Store the parent folder Xe in your root directory C:/
Then run MAIN.PY in a python terminal or IDE. Must be Python 2.7


Instructions:

AlphaViewer - Use this to view the data that has been processed and to process new data

Refresh: refreshes alpha frame list
Display Data: displays selected alpha frame
Save Table to Clip: saves currently displayed table to clipboard (try pasting into excel table)
Delete Frame: delete currently selected frame

New Entry: prompts user to select waveform containing file, puts processed data into currently displayed waveform
Create New Alpha DF: if checked, puts new waveform in new dataframe with name in adjacent box
Select Multiple Runs: if checked, "New Entry" will prompt user to select folder of runs, (folder
containing waveform containing folders- will only select "SelfTrig" runs that have not been processed yet)

Display Histogram: Displays Histogram of currently selected row in table
Delete Row: Deletes currently selected row from the table

----------------------------------------------------------

RunFitViewer - Use this to plot the data from the tables

Separation: Select a separation
Bias V: Select a bias voltage
Teflon: if checked, only shows data after Teflon reflectors were added
Postbaking: if checked, only shows data after the baking incident
Refresh: after choosing previous 4 settings, click refresh to show valid dates
Date: select which data you want to view (if there are multiple consecutive days,
this will plot data from all of them regardless of which you choose)
Plot Data: plots all data matching selections and accompanying fit

----------------------------------------------------------

Environment List:

You will need this text to create the python environment to run the app. This will make sure you have
all of the dependencies. Copy paste everything below and follow instructions. This is google-able. (Using Anaconda Spyder 2.7).

# This file may be used to create an environment using:
# $ conda create --name <env> --file <this file>
# platform: win-64
_ipyw_jlab_nb_ext_conf=0.1.0=py27_0
alabaster=0.7.12=py27_0
anaconda=2019.03=py27_0
anaconda-client=1.7.2=py27_0
anaconda-navigator=1.9.7=py27_0
anaconda-project=0.8.2=py27_0
asn1crypto=0.24.0=py27_0
astroid=1.6.5=py27_0
astropy=2.0.9=py27hc997a72_0
atomicwrites=1.3.0=py27_1
attrs=19.1.0=py27_1
babel=2.6.0=py27_0
backports=1.0=py27_1
backports.functools_lru_cache=1.5=py27_1
backports.os=0.1.1=py27_0
backports.shutil_get_terminal_size=1.0.0=py27_2
backports.shutil_which=3.5.2=py27_0
backports_abc=0.5=py27h0ec6b72_0
beautifulsoup4=4.7.1=py27_1
bitarray=0.8.3=py27h0c8e037_0
bkcharts=0.2=py27h92b6de3_0
blas=1.0=mkl
bleach=3.1.0=py27_0
blosc=1.15.0=hc65f11a_0
bokeh=1.0.4=py27_0
boto=2.49.0=py27_0
bottleneck=1.2.1=py27hc997a72_1
bzip2=1.0.6=h0c8e037_5
ca-certificates=2019.1.23=0
cdecimal=2.3=py27h0c8e037_3
certifi=2019.3.9=py27_0
cffi=1.12.2=py27hcfb25f9_1
chardet=3.0.4=py27_1
click=7.0=py27_0
cloudpickle=0.8.0=py27_0
clyent=1.2.2=py27_1
colorama=0.4.1=py27_0
comtypes=1.1.7=py27_0
conda=4.6.11=py27_0
conda-build=3.17.8=py27_0
conda-env=2.6.0=1
conda-verify=3.1.1=py27_0
configparser=3.7.3=py27_1
console_shortcut=0.1.1=3
contextlib2=0.5.5=py27h42efda5_0
cryptography=2.4.2=py27hc64555f_0
curl=7.64.0=h7a46e7a_2
cycler=0.10.0=py27h59acbbf_0
cython=0.29.6=py27hc56fc5f_0
cytoolz=0.9.0.1=py27h0c8e037_1
dask=1.1.4=py27_1
dask-core=1.1.4=py27_1
decorator=4.4.0=py27_1
defusedxml=0.5.0=py27_1
distributed=1.26.0=py27_1
docutils=0.14=py27h8652d09_0
entrypoints=0.3=py27_0
enum34=1.1.6=py27_1
et_xmlfile=1.0.1=py27h1de5d23_0
fastcache=1.0.2=py27h0c8e037_2
filelock=3.0.10=py27_0
flask=1.0.2=py27_1
freetype=2.9.1=h4d385ea_1
funcsigs=1.0.2=py27h8885ae1_0
functools32=3.2.3.2=py27_1
future=0.17.1=py27_0
futures=3.2.0=py27_0
get_terminal_size=1.0.0=h38e98db_0
gevent=1.4.0=py27h0c8e037_0
glob2=0.6=py27_1
greenlet=0.4.15=py27h0c8e037_0
grin=1.2.1=py27_4
h5py=2.9.0=py27hb721d18_0
hdf5=1.10.4=h530792d_0
heapdict=1.0.0=py27_2
html5lib=1.0.1=py27_0
icc_rt=2019.0.0=h0cc432a_1
icu=58.2=h2aa20d9_1
idna=2.8=py27_0
imageio=2.5.0=py27_0
imagesize=1.1.0=py27_0
importlib_metadata=0.8=py27_0
intel-openmp=2019.3=203
ipaddress=1.0.22=py27_0
ipykernel=4.10.0=py27_0
ipython=5.8.0=py27_0
ipython_genutils=0.2.0=py27hbe997df_0
ipywidgets=7.4.2=py27_0
isort=4.3.16=py27_0
itsdangerous=1.1.0=py27_0
jdcal=1.4=py27_0
jedi=0.13.3=py27_0
jinja2=2.10=py27_0
jpeg=9b=ha175dff_2
jsonschema=3.0.1=py27_0
jupyter=1.0.0=py27_7
jupyter_client=5.2.4=py27_0
jupyter_console=5.2.0=py27_1
jupyter_core=4.4.0=py27_0
jupyterlab=0.33.11=py27_0
jupyterlab_launcher=0.11.2=py27h28b3542_0
keyring=18.0.0=py27_0
kiwisolver=1.0.1=py27hc56fc5f_0
krb5=1.16.1=hb4d044e_6
lazy-object-proxy=1.3.1=py27h0c8e037_2
libarchive=3.3.3=h96cdc4e_0
libcurl=7.64.0=h7a46e7a_2
libiconv=1.15=hda2e4ec_7
libpng=1.6.36=h7a46e7a_0
libsodium=1.0.16=h8b3e59e_0
libssh2=1.8.0=h77a7533_4
libtiff=4.0.10=h1c3b264_2
libxml2=2.9.9=h325896a_0
libxslt=1.1.33=h803002f_0
linecache2=1.0.0=py27_0
llvmlite=0.28.0=py27hc56fc5f_0
locket=0.2.0=py27h1ca288a_1
lxml=4.3.2=py27h31b8cb8_0
lz4-c=1.8.1.2=h3cc03e0_0
lzo=2.10=h0bb7fe3_2
m2w64-gcc-libgfortran=5.3.0=6
m2w64-gcc-libs=5.3.0=7
m2w64-gcc-libs-core=5.3.0=7
m2w64-gmp=6.1.0=2
m2w64-libwinpthread-git=5.0.0.4634.697f757=2
markupsafe=1.1.1=py27h0c8e037_0
matplotlib=2.2.3=py27h263d877_0
mccabe=0.6.1=py27_1
menuinst=1.4.16=py27h0c8e037_0
mistune=0.8.4=py27h0c8e037_0
mkl=2019.3=203
mkl-service=1.1.2=py27h0b88c2a_5
mkl_fft=1.0.10=py27h44c1dab_0
more-itertools=5.0.0=py27_0
mpmath=1.1.0=py27_0
msgpack-python=0.6.1=py27hdc96acc_1
msys2-conda-epoch=20160418=1
multipledispatch=0.6.0=py27_0
navigator-updater=0.2.1=py27_0
nbconvert=5.4.1=py27_3
nbformat=4.4.0=py27hf49b375_0
networkx=2.2=py27_1
nltk=3.4=py27_1
nose=1.3.7=py27_2
notebook=5.7.8=py27_0
numba=0.43.1=py27h39f3610_0
numexpr=2.6.9=py27haac76bc_0
numpy=1.16.2=py27h5fc8d92_0
numpy-base=1.16.2=py27hb1d0314_0
numpydoc=0.8.0=py27_0
olefile=0.46=py27_0
openpyxl=2.6.1=py27_1
openssl=1.0.2r=h0c8e037_0
packaging=19.0=py27_0
pandas=0.24.2=py27hc56fc5f_0
pandoc=2.2.3.2=0
pandocfilters=1.4.2=py27_1
parso=0.3.4=py27_0
partd=0.3.10=py27_1
path.py=11.5.0=py27_0
pathlib2=2.3.3=py27_0
patsy=0.5.1=py27_0
pep8=1.7.1=py27_0
pickleshare=0.7.5=py27_0
pillow=5.4.1=py27h5b88493_0
pip=19.0.3=py27_0
pkginfo=1.5.0.1=py27_0
pluggy=0.9.0=py27_0
ply=3.11=py27_0
powershell_shortcut=0.0.1=2
prometheus_client=0.6.0=py27_0
prompt_toolkit=1.0.15=py27h3a8ec6a_0
psutil=5.5.0=py27h0c8e037_0
py=1.8.0=py27_0
pycodestyle=2.5.0=py27_0
pycosat=0.6.3=py27h0c8e037_0
pycparser=2.19=py27_0
pycrypto=2.6.1=py27h0c8e037_9
pycurl=7.43.0.2=py27hc64555f_0
pyflakes=2.1.1=py27_0
pygments=2.3.1=py27_0
pylint=1.9.2=py27_0
pyodbc=4.0.26=py27hc56fc5f_0
pyopenssl=19.0.0=py27_0
pyparsing=2.3.1=py27_0
pyqt=5.6.0=py27h6e61f57_6
pyreadline=2.1=py27_1
pyrsistent=0.14.11=py27h0c8e037_0
pysocks=1.6.8=py27_0
pytables=3.5.1=py27h6a9b274_0
pytest=4.3.1=py27_0
python=2.7.16=hcb6e200_0
python-dateutil=2.8.0=py27_0
python-libarchive-c=2.8=py27_6
pytz=2018.9=py27_0
pywavelets=1.0.2=py27hc997a72_0
pywin32=223=py27h0c8e037_1
pywinpty=0.5.5=py27_1000
pyyaml=5.1=py27h0c8e037_0
pyzmq=18.0.0=py27hc56fc5f_0
qt=5.6.2=vc9hc26998b_12
qtawesome=0.5.7=py27_1
qtconsole=4.4.3=py27_0
qtpy=1.7.0=py27_1
requests=2.21.0=py27_0
rope=0.12.0=py27_0
ruamel_yaml=0.15.46=py27h0c8e037_0
scandir=1.10.0=py27h0c8e037_0
scikit-image=0.14.2=py27hc56fc5f_0
scikit-learn=0.20.3=py27hf381715_0
scipy=1.2.1=py27h4c3ab11_0
seaborn=0.9.0=py27_0
send2trash=1.5.0=py27_0
setuptools=40.8.0=py27_0
simplegeneric=0.8.1=py27_2
singledispatch=3.4.0.3=py27h3f9d112_0
sip=4.18.1=py27hc56fc5f_2
six=1.12.0=py27_0
snappy=1.1.7=he46498f_3
snowballstemmer=1.2.1=py27h28d3bf7_0
sortedcollections=1.1.2=py27_0
sortedcontainers=2.1.0=py27_0
soupsieve=1.8=py27_0
sphinx=1.8.5=py27_0
sphinxcontrib=1.0=py27_1
sphinxcontrib-websupport=1.1.0=py27_1
spyder=3.3.3=py27_0
spyder-kernels=0.4.2=py27_0
sqlalchemy=1.3.1=py27h0c8e037_0
sqlite=3.27.2=h0c8e037_0
ssl_match_hostname=3.7.0.1=py27_0
statsmodels=0.9.0=py27hc997a72_0
subprocess32=3.5.3=py27h0c8e037_0
sympy=1.3=py27_0
tblib=1.3.2=py27h8ae915c_0
terminado=0.8.1=py27_1
testpath=0.4.2=py27_0
tk=8.6.8=h0c8e037_0
toolz=0.9.0=py27_0
tornado=5.1.1=py27h0c8e037_0
tqdm=4.31.1=py27_1
traceback2=1.4.0=py27_0
traitlets=4.3.2=py27h1b1b3a5_0
typing=3.6.6=py27_0
unicodecsv=0.14.1=py27h0bf7bb0_0
unittest2=1.1.0=py27_0
urllib3=1.24.1=py27_0
vc=9=h7299396_1
vs2008_runtime=9.00.30729.1=hfaea7d5_1
vs2015_runtime=14.15.26706=h3a45250_0
wcwidth=0.1.7=py27hb1a0d82_0
webencodings=0.5.1=py27_1
werkzeug=0.14.1=py27_0
wheel=0.33.1=py27_0
widgetsnbextension=3.4.2=py27_0
win_inet_pton=1.1.0=py27_0
win_unicode_console=0.5=py27hc037021_0
wincertstore=0.2=py27hf04cefb_0
winpty=0.4.3=4
wrapt=1.11.1=py27h0c8e037_0
xlrd=1.2.0=py27_0
xlsxwriter=1.1.5=py27_0
xlwings=0.15.4=py27_0
xlwt=1.3.0=py27h2271735_0
xz=5.2.4=h3cc03e0_4
yaml=0.1.7=h3e6d941_2
zeromq=4.3.1=h2880e7c_3
zict=0.1.4=py27_0
zipp=0.3.3=py27_1
zlib=1.2.11=h3cc03e0_3
zstd=1.3.7=h1b0e4d7_0
