import math
import os
import random
import argparse
from datetime import datetime

import numpy as np
import tensorflow as tf
from tensorflow import keras

from utils.structure import getModel
from utils.dataset import gen, get_dataset_from_csv
from utils.record import saveConfig, drawFig


print(tf.test.gpu_device_name())
AUTOTUNE = tf.data.experimental.AUTOTUNE


def main(args):
    # load dataset
    random.seed(args.seed)
    np.random.seed(args.seed)
    tf.random.set_seed(args.seed)

    print('---------------loading dataset---------------')
    input_size = [args.alength, args.alength, args.alength, args.ndims]
    train_set = get_dataset_from_csv(args.train_set, args.datapath)
    test_set = get_dataset_from_csv(args.test_set, args.datapath)

    print('dataset loaded, train set size', len(train_set), 'test set size', len(test_set))

    print('---------------training---------------')

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    savingPath = os.path.join(args.savingPath, timestamp)
    print("output path:", savingPath)

    if not os.path.exists(savingPath):
        os.makedirs(savingPath)

    saveConfig(savingPath, args)

    examine = 'accuracy'  # binary_accuracy accuracy
    monitor = f'val_{examine}'

    earlystopper = tf.keras.callbacks.EarlyStopping(
        monitor=monitor, patience=20, verbose=1)

    save_best = tf.keras.callbacks.ModelCheckpoint(
    filepath=os.path.join(savingPath, "best_model.weights.h5"), 
    monitor=monitor, 
    verbose=1, 
    save_best_only=False, 
    mode='max' if examine == 'accuracy' else 'min'  # Ensure mode matches the metric
)    
    save_latest = tf.keras.callbacks.ModelCheckpoint(
    filepath=os.path.join(savingPath, "latest_model.weights.h5"), 
    save_freq='epoch',  # Save every epoch
    verbose=1
    )

    print(f"Saving models to: {savingPath}")


    model = getModel(args.model, input_size)

    print("Training completed. Models should be saved.")

    if args.weights:
        print("Using trained weights from ")
        model.load_weights(args.weights)

    # model.summary()
    lossfun = keras.losses.CategoricalCrossentropy(from_logits=False)

    model.compile(optimizer=keras.optimizers.Adam(learning_rate=args.lr),
                loss=lossfun,
                metrics=[examine, tf.keras.metrics.AUC()])

    history = model.fit(
        gen(train_set, args.batch, shuffle=True, augment=args.augment),
        epochs=args.epoch,
        steps_per_epoch=math.ceil(len(train_set)/args.batch),
        validation_data=gen(test_set, args.batch, shuffle=False, augment=False),
        validation_steps=math.ceil(len(test_set)/args.batch),
        callbacks=[save_best, earlystopper, save_latest],
        verbose=1
    )
    drawFig(history, examine, monitor, savingPath)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', choices=['Resnet3D', 'DenseNet3D'], default='DenseNet3D', help='Backbone of the network')
    parser.add_argument('--datapath', type=str, default='./data/example_dataset/distance', help='Path to tensors')
    parser.add_argument('--weights', type=str, default='', help='Continue training based on weights')
    parser.add_argument('--savingPath', type=str, default='./models/example_model', help='Path to save models')

    parser.add_argument('--train_set', type=str, default='./data/example_dataset/part_0_train.csv', help='Path to train set csv file')
    parser.add_argument('--test_set', type=str, default='./data/example_dataset/part_0_val.csv', help='Path to val set csv file')

    parser.add_argument('--augment', type=bool, default=True, help='Use data augmentation')
    parser.add_argument('--batch', default=32, type=int)
    parser.add_argument('--alength', default=64, type=int, help='Tensor side length')
    parser.add_argument('--ndims', default=8, type=int, help='Tensor dims')
    parser.add_argument('--seed', default=2032, type=int, help='Random seeds')
    parser.add_argument('--epoch', default=200, type=int, help='Epoch to train')
    parser.add_argument('--lr', default=1e-4, type=float, help='Init learning rate')
    args, unknown = parser.parse_known_args()
    main(args)
