# GP-Network-Classification
This repository is for work on developing a framework for classification with network inputs, using Gaussian process models. Suggested organizational structure is to have subdirectories for

- data (but do NOT store either proprietary or large data sets)
- docs (e.g., LaTeX source files, with accompanying style files, etc.)
- code
- results (e.g., key figures, etc. Due to LaTeX conventions, you might have a separte img dir within docs)

Some notes on best-practice usage:

- Git works best for text files. It stores files and then updates then. It does this easiest with text.
- If, for example, you constantly update an image file, it will essentially store multiple copies of that file, which is inefficient. (Note: You can use a "gitignore" command on such files to suppress this.) You make changes locally to a file, then push that file to the git repo and commit them. How often you do this is a matter of choice, and driven to some extent by context. With software development, people often commit even fairly small changes, once the code is established. With text, you might commit larger changes. One key thing that git allows you to do is to use a "diff" command to see what a collaborator has changed (or you yourself, if you've forgotten!). So think about what order of magnitude changes you'd like to see in the diff command.
- In general, if you can recreate a file locally (i.e., on your computer), do NOT store it on git. Unless, of course, you judge it to be simply too expensive to want to recreate. So, for example, .tex files are ok, but the output (i.e., .aux, .log, .dvi, .pdf, etc.) are unnecessary.
