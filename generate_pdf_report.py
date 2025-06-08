# generate_pdf_report.py
# Converts an HTML file to a PDF file using WeasyPrint.
#
# Prerequisites:
#   Install WeasyPrint: pip install WeasyPrint
#   WeasyPrint may also have system-level dependencies for font rendering (e.g., Pango, Cairo, GDK-PixBuf).
#   Refer to WeasyPrint documentation for installation details: https://weasyprint.readthedocs.io/en/stable/install.html
#
# Usage:
#   python generate_pdf_report.py <input_html_file> <output_pdf_file>
#
# Example:
#   python generate_pdf_report.py comprehensive_research_summary.html comprehensive_research_summary.pdf

import argparse
import os
from weasyprint import HTML, CSS

def convert_html_to_pdf(html_filepath, pdf_filepath):
    """
    Converts the given HTML file to a PDF file.

    Args:
        html_filepath (str): The path to the input HTML file.
        pdf_filepath (str): The path where the output PDF file will be saved.
    """
    try:
        print(f"Starting conversion of '{html_filepath}' to '{pdf_filepath}'...")

        # Ensure the input HTML file exists
        if not os.path.exists(html_filepath):
            print(f"Error: Input HTML file not found at '{html_filepath}'")
            return

        # Base URL for resolving relative paths (e.g., for images, CSS)
        base_url = os.path.dirname(os.path.abspath(html_filepath))

        # Basic CSS for page size, margins, and default font
        css = CSS(string='''
            @page { size: A4; margin: 1in; }
            body { font-family: sans-serif; }
        ''')

        # HTML object
        html_doc = HTML(filename=html_filepath, base_url=base_url)
        
        # Write PDF with the defined CSS
        html_doc.write_pdf(pdf_filepath, stylesheets=[css])

        print(f"Successfully converted '{html_filepath}' to '{pdf_filepath}'")

    except Exception as e:
        print(f"An error occurred during PDF conversion: {e}")
        print("Please ensure WeasyPrint and its dependencies are correctly installed.")
        print("See WeasyPrint documentation: https://weasyprint.readthedocs.io/en/stable/install.html")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert an HTML file to PDF using WeasyPrint.")
    parser.add_argument("input_html_file", help="Path to the input HTML file.")
    parser.add_argument("output_pdf_file", help="Path to save the output PDF file.")

    args = parser.parse_args()

    convert_html_to_pdf(args.input_html_file, args.output_pdf_file)
