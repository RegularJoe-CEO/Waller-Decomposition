FROM python:3.12
COPY . /app
RUN pip install -e .
CMD ["streamlit", "run", "ui/app.py"]