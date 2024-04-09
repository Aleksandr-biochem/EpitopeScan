from setuptools import setup

setup(
    name='EpitopeScan',
    version='0.1.0',    
    description='Python3 toolset for the analysis of SARS-CoV-2 genomes reporting mutations in immunogenic epitopes of interest',
    url='https://github.com/Aleksandr-biochem/EpitopeScan',
    author='Aleksandr Kovalenko, Sebastien Viatte',
    author_email='alekskov1102@gmail.com',
    license='MIT Licenses',
    packages= ['EpitopeScan', 'EpitopeScan.utils'],
    package_data={'EpitopeScan': ['reference_sequences/*']},
    install_requires=['numpy>=1.24.0',
                      'blosum>=2.0.2',
                      'pandas>=2.0.1',
                      'streamlit>=1.22.0',
                      'plotly>=5.14.1'],
    entry_points={
        'console_scripts': [
            'EpitopeScan = EpitopeScan:EpitopeScan.main',
            'EpitopeScanGUI = EpitopeScan:EpitopeScanGUI.run_gui'
        ]
    },
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research'
    ],
)