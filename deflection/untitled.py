import click
import yaml


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
# @click.option('--n_cpu', default=1, help='Number of CPUs')
# @click.argument('a', default=0)



def main(cfg):
    with open(cfg, 'r') as stream:
        cfg = yaml.full_load(stream)
        
    print(stream.name)
    s = stream.name # "configs/test.yaml"
    
    config_name = s[len(s) - s[::-1].find('/'):s.find('.yaml')]
    print(config_name)




        
        
    

if __name__ == '__main__':
    
    
    main()