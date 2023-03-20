import sys
import datetime
__version__ = "1.1.0.0"
__release_yearMonthDate__ = datetime.date.today()


def main():
    import pyinstaller_versionfile
    pyinstaller_versionfile.create_versionfile(
        output_file="version.rc",
        version=__version__,
        company_name="Vladislav Powerman and the Institute of Earths Crust, SB RAS",
        file_description="Dezirteer",
        internal_name="Dezirteer",
        legal_copyright="©Vladislav Powerman and the Institute of Earths Crust. All rights reserved.",
        original_filename="dezirteer.exe",
        product_name="Dezirteer"
    )
if __name__ == '__main__':
    main()
