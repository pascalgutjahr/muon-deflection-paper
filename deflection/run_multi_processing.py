import os 
import click

@click.command()
@click.argument('folder')

def main(folder):
    '''Run multi_processing script with several configs.
    Parameters
    ----------
        folder: string
            folder to settings, example: configs/settings
    '''
    
    print(folder)
    
    for file in os.listdir(folder):
        os.system('python multi_processing.py {}/{}'.format(folder, file))
        
        
if __name__ == '__main__':
    
    main()