import { Observable } from 'rxjs';

export interface DatasetInfoServiceInterface {
  datasetInfo: Observable<string[]>;

  setDatasetInfo(): Promise<void>;
}