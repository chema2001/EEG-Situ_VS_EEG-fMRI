# Documentação do Brainflow: https://brainflow.readthedocs.io/en/stable/UserAPI.html

'''
Este scrip é executado a partir do terminal de forma a definir os elétrodos onde se quer adquirir o EEG

Para correr o script é preciso primeiro definir o diretório onde ele está:
- cd path

De seguida:
- python acquisition_brainflow.py --eeg cyton --serial-port COM6 --channels FP1,FP2,...

Para testar o script sem nenhuma placa conectada:
- python acquisition_brainflow.py --eeg test
'''


import keyboard

import argparse
import logging
import numpy as np

from brainflow.board_shim import BoardShim, BrainFlowInputParams, BoardIds
from brainflow.data_filter import DataFilter, FilterTypes, DetrendOperations

def eeg_systemId(_id):
    
    # Seleção da placa de aquisição a partir do ID da mesma na API do brainflow.
    # Com este ID, a API faz o setup com as predefinições necessárias para adquirir sinal da placa

    '''
    - A placa daisy é uma extensão da cyton que permite adicionar mais 8 elétrodos. 
    - Em princípio isto não irá ser preciso a não ser que queiramos ter melhor resolução espacial.
    - Filtros espaciais (ICA, CSP) por vezes usados no pré-processamento (retirar ruído) e na análise de ERPs, 
      funcionam melhor com um maior número de elétrodos por isso não é má ideia considerarmos isto 

    - A placa sintética serve para fazer testes e usa sinais simulados na biblioteca
    '''
    if(_id == "cyton"):
        board_id = BoardIds.CYTON_BOARD
        isCyton = True 
    elif(_id == "daisy"):
        board_id = BoardIds.CYTON_DAISY_BOARD 
        isCyton = True
    elif(_id == "test"):
        board_id = BoardIds.SYNTHETIC_BOARD
        isCyton = False
    return board_id, isCyton

def cyton_channels_sep(channel_string):

    # Separação da string de elétrodos nas strings individuais correspondentes a cada elétrodo 

    individual_chs = channel_string.split(',')
    return individual_chs

def main():
    BoardShim.enable_dev_board_logger()
    logging.basicConfig(level=logging.DEBUG) # Ativa as mensagens log do brainflow para fazer debug

    parser = argparse.ArgumentParser() # Permite definir os parâmetros através da consola sem ser necessário ter o editor aberto
    parser.add_argument('--eeg', type=str, # Este é o único parâmetro necessário colocar para testar (SYNTHETIC_BOARD)
                        help='EEG headband, check docs to get a list of supported boards: https://brainflow.readthedocs.io/en/stable/SupportedBoards.html#supported-boards-label',
                        required=True, default='')
    parser.add_argument('--channels', type=cyton_channels_sep,
                        help='Headbad channels that will be used. For example, write in order and seperated by a comma from pins NP1 to NP8: [Fp1, Fp2,...]',
                        required=False, default='')
    parser.add_argument('--serial-port', type=str, help='Serial port', required=False, default='')
    args = parser.parse_args()

    params = BrainFlowInputParams()
    params.serial_port = args.serial_port 
    
    '''
    O serial port (no caso da cyton e daisy) é o único parâmetro necessário colocar diretamente como argumento
    para inicializar a board com as predefinições do brainflow. Os outros dois são usados à parte 
    '''
    try:  
        board_id, isCyton = eeg_systemId(args.eeg)
        if (isCyton):
            cyton_ch = args.channels
        else:
            cyton_ch = ""
        print("\nBoard description: \n")
        print(BoardShim.get_board_descr(board_id,0))
        print("\n")
        board = BoardShim(board_id, params) # Inicialização da board
        board.prepare_session()
        board.start_stream(450000) # O parâmetro corresponde ao tamanho do Buffer. O valor é o default do brainflow
        eeg_channels = BoardShim.get_eeg_channels(board_id) # O brainflow adquire dados de 24 colunas (no caso da cyton) mas apenas 8 delas correspondem aos dados do EEG
        sampling_rate = BoardShim.get_sampling_rate(board_id)

        if(len(cyton_ch) > 1): 
            # Caso os nomes dos elétrodos sejam inseridos como parâmetro, os mesmos são guardados numa variável
            channel_labels = cyton_ch
            eeg_channels = eeg_channels[0:len(channel_labels)]
        else:
            # Caso não sejam inseridos, os nomes dos elétrodos predifinidos no brainflow para a placa escolhida serão usados para não dar erro
            channel_labels = BoardShim.get_eeg_names(board_id)

        # Variável com os dados da placa. 
        # O parâmetro define a quantidade de amostras a retirar do buffer. Neste caso retira as amostras de 1s 
        while(True):
            data = board.get_current_board_data(sampling_rate)
            print("\n")
            print(data[eeg_channels]) 
            if keyboard.is_pressed('0'): # Clicar no 0 para parar o stream.
                stopStream(board)
                break


    except BaseException:
        print("--------------------------------------------------------------------------------------------------------------------------")
        logging.warning('Exception', exc_info=True)
        print("--------------------------------------------------------------------------------------------------------------------------")
        
    finally:
        logging.info('End')
        if board.is_prepared():
            stopStream(board)

def stopStream(board): # Termina a sessão corretamente
    logging.info('Releasing session')
    board.release_session()

if __name__ == '__main__':
    main()
