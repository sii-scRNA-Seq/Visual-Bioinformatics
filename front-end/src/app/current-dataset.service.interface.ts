import { Observable } from 'rxjs';

import { DatasetInfo } from './dataset-info';

export interface CurrentDatasetServiceInterface {

  currentDataset: Observable<DatasetInfo>;

  setCurrentDataset(newDataset: string): void;
}
