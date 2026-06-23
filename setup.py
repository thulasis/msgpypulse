from setuptools import setup, find_packages

setup(
    name="msgpypulse",
    version="4.1.1",
    packages=find_packages(),  # This will find the msgpypulse/ directory
    install_requires=[
        
    ],
    entry_points={
        'console_scripts': [
            'msgpypulse=scripts.msgpypulse_cli:main',        # Calls main() in msgpypulse_cli.py
            'msgpypulse-gui=scripts.msgpypulse_gui:main', # Calls main() in msgpypulse_gui.py
        ],
    },
    author="Tulasi Relangi & Harrison Hall, 2025",
    description="MSGF+ downstream analysis module",
)