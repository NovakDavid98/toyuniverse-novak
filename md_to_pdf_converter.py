import markdown2
from weasyprint import HTML, CSS
import os

MARKDOWN_FILE = "simulation.md"
HTML_FILE = "temp_simulation.html" # Temporary HTML file
PDF_FILE = "simulation.pdf"

CSS_STYLE = """
@import url('https://fonts.googleapis.com/css2?family=Noto+Sans:ital,wght@0,400;0,700;1,400&family=Noto+Serif:ital,wght@0,400;0,700;1,400&family=Noto+Mono&display=swap');

body {
    font-family: 'Noto Serif', serif;
    line-height: 1.6;
    color: #333333;
    margin: 2cm;
    background-color: #FCFCFC; /* Very light off-white background */
}

h1, h2, h3, h4, h5, h6 {
    font-family: 'Noto Sans', sans-serif;
    color: #005A9C; /* Dark professional blue */
    margin-top: 1.8em;
    margin-bottom: 0.6em;
    line-height: 1.3;
}

h1 {
    font-size: 2.4em;
    text-align: center;
    border-bottom: 2px solid #005A9C;
    padding-bottom: 0.3em;
    margin-bottom: 1.2em;
    margin-top: 0;
}

h2 {
    font-size: 1.9em;
    border-bottom: 1px solid #007ACC; /* Lighter accent blue */
    padding-bottom: 0.2em;
}

h3 {
    font-size: 1.5em;
    color: #007ACC;
}

h4 {
    font-size: 1.2em;
    color: #1E88E5; /* A slightly brighter blue for H4 */
}

p {
    margin-bottom: 1em;
    text-align: left; /* Justify can sometimes look odd in PDFs from HTML */
}

a {
    color: #007ACC;
    text-decoration: none;
}

a:hover {
    text-decoration: underline;
    color: #005A9C;
}

strong, b {
    font-weight: bold;
    color: #222222; /* Slightly darker for emphasis */
}

em, i {
    font-style: italic;
    color: #444444;
}

ul, ol {
    margin-bottom: 1em;
    padding-left: 1.8em;
}

li {
    margin-bottom: 0.3em;
}

pre {
    background-color: #F0F4F8; /* Light grayish blue */
    border: 1px solid #D0DAE0;
    border-left: 5px solid #007ACC;
    padding: 1em;
    border-radius: 5px;
    overflow-x: auto;
    font-family: 'Noto Mono', monospace;
    font-size: 0.9em;
    line-height: 1.4;
    margin-top: 1em;
    margin-bottom: 1.5em;
    white-space: pre-wrap; /* Allow wrapping for long lines */
    word-wrap: break-word; /* Break words if necessary */
}

code {
    font-family: 'Noto Mono', monospace;
    background-color: #E0E8F0; /* Slightly darker for inline */
    padding: 0.2em 0.4em;
    border-radius: 3px;
    font-size: 0.85em;
    color: #2C3E50; /* Dark blue-gray for inline code text */
}

pre code { /* Reset for code inside pre */
    background-color: transparent;
    padding: 0;
    border-radius: 0;
    font-size: inherit; /* Inherit from pre */
    color: inherit;
}

blockquote {
    border-left: 4px solid #007ACC;
    padding-left: 1.5em;
    margin-left: 0;
    margin-right: 0;
    margin-top: 1.5em;
    margin-bottom: 1.5em;
    color: #555555;
    font-style: italic;
    background-color: #F8F9FA; /* Very light background for blockquotes */
}

blockquote p {
    margin-bottom: 0.5em;
}

img {
    max-width: 90%; /* Avoid images touching the very edge of content area */
    height: auto;
    display: block;
    margin-left: auto;
    margin-right: auto;
    margin-top: 1.5em;
    margin-bottom: 1.5em;
    border: 1px solid #DDDDDD;
    border-radius: 4px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.05);
}

hr {
    border: 0;
    height: 1px;
    background: #CCCCCC;
    margin: 2em 0;
}

/* Table styling */
table {
    width: 100%;
    border-collapse: collapse;
    margin-bottom: 1.5em;
    font-size: 0.9em;
}

th, td {
    border: 1px solid #DDDDDD;
    padding: 0.75em;
    text-align: left;
}

th {
    background-color: #E9ECEF; /* Light gray for table headers */
    font-family: 'Noto Sans', sans-serif;
    font-weight: bold;
    color: #005A9C;
}

/* Page numbers */
@page {
    @bottom-right {
        content: "Page " counter(page) " of " counter(pages);
        font-family: 'Noto Sans', sans-serif;
        font-size: 0.8em;
        color: #777777;
    }
}
"""

def create_pdf():
    """
    Converts the markdown file to a styled PDF.
    """
    try:
        # Read markdown content
        with open(MARKDOWN_FILE, 'r', encoding='utf-8') as f:
            md_content = f.read()

        # Convert markdown to HTML
        # Using 'fenced-code-blocks' and 'tables' extras for common markdown features
        html_content = markdown2.markdown(md_content, extras=['fenced-code-blocks', 'tables', 'footnotes', 'smarty-pants'])

        # Full HTML document structure
        full_html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <title>{os.path.splitext(MARKDOWN_FILE)[0]}</title>
        </head>
        <body>
            {html_content}
        </body>
        </html>
        """
        
        css = CSS(string=CSS_STYLE)
        
        # Base URL is important for resolving relative paths for images
        base_url = os.path.dirname(os.path.abspath(MARKDOWN_FILE))

        # Create HTML object
        html_doc = HTML(string=full_html, base_url=base_url)
        
        # Write PDF
        html_doc.write_pdf(PDF_FILE, stylesheets=[css])
        
        print(f"Successfully converted '{MARKDOWN_FILE}' to '{PDF_FILE}'")

    except FileNotFoundError:
        print(f"Error: Markdown file '{MARKDOWN_FILE}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    create_pdf()
