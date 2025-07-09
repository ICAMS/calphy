def test_ex07():
    import calphy
    from calphy.postprocessing import read_report

# def test_example07_first_cell_runs():
#     notebook_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../examples/example_07/analysis.ipynb'))

#     with open(notebook_path, "r", encoding="utf-8") as f:
#         nb = nbformat.read(f, as_version=4)

#     # Filter just the first code cell
#     first_code_cell = next((cell for cell in nb.cells if cell.cell_type == "code"), None)

#     if not first_code_cell:
#         pytest.fail("No code cell found in example07.ipynb")

#     # Create a minimal notebook with only the first cell
#     nb_single = nbformat.v4.new_notebook()
#     nb_single.cells = [first_code_cell]

#     try:
#         client = NotebookClient(nb_single, timeout=60, kernel_name="pyiron")
#         client.execute()
#     except Exception as e:
#         pytest.fail(f"First code cell in example07.analysis.ipynb failed: {e}")
